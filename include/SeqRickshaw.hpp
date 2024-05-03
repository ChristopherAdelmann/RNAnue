#ifndef SEQRICKSHAW_H
#define SEQRICKSHAW_H

#include "Logger.hpp"
#include "Utility.hpp"

// openMP
#include <omp.h>

#include <algorithm>
#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <ranges>
#include <set>
#include <thread>
#include <tuple>
#include <vector>

// boost
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan
#include <numeric>
#include <ranges>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/slice.hpp>

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

typedef std::tuple<std::string, int, std::pair<int, int>, std::size_t> State;
typedef std::vector<State> States;

typedef std::map<std::pair<int, char>, std::tuple<int, int, int, int>> LookupTable;

typedef std::map<std::pair<std::string, std::string>, LookupTable> Adapters;

typedef std::vector<fs::path> PathVector;
typedef std::pair<PathVector, PathVector> PathVectorPair;

// typedef std::tuple<seqan3::field::id,seqan3::field::seq,seqan3::field:qual>

// typedef std::map<std::pair<int,string::string>
//

using seqan3::operator""_dna5;
using seqan3::operator""_dna4;
class SeqRickshaw {
   private:
    const std::string readtype;

    const bool trimPolyG;

    const std::string adpt5f;
    const std::string adpt5r;
    const std::string adpt3f;
    const std::string adpt3r;
    const double missMatchRateTrim;
    const int minOverlapTrim;

    const int minMeanPhread;
    const int minLen;

    const int minWindowPhread;
    const int windowTrimSize;

    const int minOverlapMerge;
    const double missMatchRateMerge;

    const size_t threads;
    const size_t chunkSize;

    struct SingleEndFastqChunk {
        std::vector<seqan3::sequence_file_input<>::record_type> records;
    };

    struct PairedEndFastqChunk {
        std::vector<seqan3::sequence_file_input<>::record_type> recordsFwd;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsRev;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsMergedRes;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsSnglFwdRes;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsSnglRevRes;
    };

    struct TrimConfig {
       public:
        enum Mode { FIVE_PRIME, THREE_PRIME };

        static seqan3::align_cfg::method_global alignmentConfigFor(Mode mode);
    };

    struct Adapter {
        const seqan3::dna5_vector sequence;
        const double maxMissMatchFraction;
        const TrimConfig::Mode trimmingMode;
    };

    const std::vector<Adapter> loadAdapters(const std::string &filenameOrSequence,
                                            const TrimConfig::Mode trimmingMode);

    const bool passesFilters(const auto &record);

    template <typename record_type>
    void trimWindowedQuality(record_type &record);

    template <typename record_type>
    void trimAdapter(const Adapter &adapter, record_type &record);

    template <typename record_type>
    void trim3PolyG(record_type &record);

    template <typename record_type>
    const std::optional<record_type> mergeRecordPair(const record_type &record1,
                                                     const record_type &record2);

    template <typename record_type, typename result_type>
    const record_type constructMergedRecord(
        const record_type &record1, const record_type &record2,
        const seqan3::alignment_result<result_type> &alignmentResult);

    void processSingleEnd(pt::ptree sample);
    void processPairedEnd(pt::ptree sample);

    void processSingleEndRecordChunk(SeqRickshaw::SingleEndFastqChunk &chunk,
                                     const std::vector<Adapter> &adapters5,
                                     const std::vector<Adapter> &adapters3);
    void processSingleEndFileInChunks(std::string const &recInPath, std::string recOutPath,
                                      const std::vector<Adapter> &adapters5,
                                      const std::vector<Adapter> &adapters3, size_t chunkSize,
                                      size_t numThreads);

    void processPairedEndRecordChunk(PairedEndFastqChunk &chunk,
                                     const std::vector<Adapter> &adapters5f,
                                     const std::vector<Adapter> &adapters3f,
                                     const std::vector<Adapter> &adapters5r,
                                     const std::vector<Adapter> &adapters3r);
    void processPairedEndFileInChunks(
        std::string const &recFwdInPath, std::string const &recRevInPath,
        std::string const &mergedOutPath, std::string const &snglFwdOutPath,
        std::string const &snglRevOutPath, const std::vector<Adapter> &adapters5f,
        const std::vector<Adapter> &adapters3f, const std::vector<Adapter> &adapters5r,
        const std::vector<Adapter> &adapters3r, size_t chunkSize, size_t numThreads);

   public:
    SeqRickshaw(const po::variables_map &params);

    void start(pt::ptree sample);
};
#endif
