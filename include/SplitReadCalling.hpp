#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

#include <mutex>

// openMP
#include <omp.h>

#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include <bitset>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <optional>

// seqan3
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

#include <seqan3/io/sam_file/all.hpp>

#include <seqan3/io/sequence_file/input.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

// filters
#include "ScoringMatrix.hpp"
#include "Traceback.hpp"
#include "Logger.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

using seqan3::operator""_tag;
using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;
using seqan3::operator""_dna15;
using seqan3::operator""_dna4;
using seqan3::get;

// overload struct to
template <>
struct seqan3::sam_tag_type<"XX"_tag>
{
    using type = int32_t;
};
template <>
struct seqan3::sam_tag_type<"XY"_tag>
{
    using type = int32_t;
};
template <>
struct seqan3::sam_tag_type<"XJ"_tag>
{
    using type = int32_t;
};
template <>
struct seqan3::sam_tag_type<"XH"_tag>
{
    using type = int32_t;
};

template <>
struct seqan3::sam_tag_type<"XM"_tag>
{
    using type = int32_t;
}; // matches in alignment
template <>
struct seqan3::sam_tag_type<"XL"_tag>
{
    using type = int32_t;
}; // length of alignment
template <>
struct seqan3::sam_tag_type<"XN"_tag>
{
    using type = int32_t;
};

template <>
struct seqan3::sam_tag_type<"XS"_tag>
{
    using type = std::string;
};

typedef std::pair<uint32_t, uint32_t> ReadPos;
typedef std::pair<uint64_t, uint64_t> GenomePos;
typedef std::vector<seqan3::cigar> CigarSplt;

// introduce record_type
using types = seqan3::type_list<
    std::string,
    seqan3::sam_flag,
    std::optional<int32_t>,
    std::optional<int32_t>,
    std::optional<uint8_t>,
    std::vector<seqan3::cigar>,
    std::span<seqan3::dna5>,
    seqan3::sam_tag_dictionary>;

using types_as_ids = seqan3::fields<
    seqan3::field::id,
    seqan3::field::flag,
    seqan3::field::ref_id,
    seqan3::field::ref_offset,
    seqan3::field::mapq,
    seqan3::field::cigar,
    seqan3::field::seq,
    seqan3::field::tags>;

using SamRecord = seqan3::sam_record<types, types_as_ids>;
using ComplResult = std::tuple<int, int, double, double, std::vector<char>, std::vector<char>>;

typedef std::vector<
    std::tuple<
        std::string,
        seqan3::sam_flag,
        std::optional<int32_t>,
        std::optional<int32_t>,
        std::vector<seqan3::cigar>,
        seqan3::dna5_vector,
        seqan3::sam_tag_dictionary>>
    Splts;

typedef std::pair<size_t, size_t> NucleotidePairPositions;
typedef std::pair<NucleotidePairPositions, NucleotidePairPositions> NucleotidePositionsWindow;
typedef std::pair<seqan3::dna5_vector, seqan3::dna5_vector> NucleotideWindowPair;

static const std::map<NucleotideWindowPair, size_t> crosslinkingScoringScheme = {
    {{"TA"_dna5, "AT"_dna5}, 3}, // Preffered pyrimidine crosslinking
    {{"TG"_dna5, "GT"_dna5}, 3},
    {{"TC"_dna5, "CT"_dna5}, 3},
    {{"AT"_dna5, "TA"_dna5}, 3},
    {{"GT"_dna5, "TG"_dna5}, 3},
    {{"CT"_dna5, "TC"_dna5}, 3},
    {{"CA"_dna5, "AC"_dna5}, 2}, // Non-preffered pyrimidine crosslinking
    {{"CG"_dna5, "GC"_dna5}, 2},
    {{"AC"_dna5, "CA"_dna5}, 2},
    {{"GC"_dna5, "CG"_dna5}, 2},
    {{"TA"_dna5, "GT"_dna5}, 1}, // Wobble base pairs crosslinking
    {{"TG"_dna5, "GT"_dna5}, 1},
    {{"GT"_dna5, "TA"_dna5}, 1},
    {{"GT"_dna5, "TG"_dna5}, 1}};
typedef struct
{
    seqan3::dna5_vector forwardWindowNucleotides;
    seqan3::dna5_vector reverseWindowNucleotides;
    std::pair<size_t, size_t> forwardWindowPositions;
    std::pair<size_t, size_t> reverseWindowPositions;
    bool isInterFragment;
} InteractionWindow;

typedef struct
{
    double energy;
    double normCrosslinkingScore;
} HybridizationResult;
class SplitReadCalling
{
private:
    po::variables_map params;
    std::vector<std::tuple<std::string>> stats;

    int readscount;
    int alignedcount;
    int splitscount;
    int msplitscount;
    int nsurvivedcount;

    double findCrosslinkingSites(
        std::span<seqan3::dna5> &seq1,
        std::span<seqan3::dna5> &seq2,
        std::vector<seqan3::dot_bracket3> &dotbracket);
    std::optional<InteractionWindow> getContinuosNucleotideWindows(
        std::span<seqan3::dna5> &seq1,
        std::span<seqan3::dna5> &seq2,
        NucleotidePositionsWindow positionsPair);

public:
    SplitReadCalling();
    SplitReadCalling(po::variables_map params);

    // iterate through reads
    void iterate(std::string matched, std::string splits, std::string multsplits);
    void process(auto &splitrecords, auto &splitsfile, auto &multsplitsfile);
    void distribute(auto &subrecords, auto &splits, auto &msplits);

    void filterSegments(auto &splitrecord, std::optional<int32_t> &refOffset,
                        std::vector<seqan3::cigar> &cigar, std::span<seqan3::dna5> &seq,
                        seqan3::sam_tag_dictionary &tags, std::vector<SamRecord> &curated);

    void addFilterToSamRecord(SamRecord &rec, std::pair<float, float> filters);
    void addComplementarityToSamRecord(SamRecord &rec1, SamRecord &rec2, TracebackResult &res);
    void addHybEnergyToSamRecord(SamRecord &rec1, SamRecord &rec2, HybridizationResult &hyb);

    void writeSamFile(auto &samfile, std::vector<std::pair<SamRecord, SamRecord>> &splits);

    TracebackResult complementarity(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2);
    std::optional<HybridizationResult> hybridize(std::span<seqan3::dna5> &seq1, std::span<seqan3::dna5> &seq2);
    void createDir(fs::path path);

    void progress(std::ostream &out);

    int countSamEntries(std::string file, std::string command);
    std::vector<std::vector<fs::path>> splitInputFile(std::string matched, std::string splits, int entries);

    std::string addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keys);
    void start(pt::ptree sample);
};
