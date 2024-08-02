#pragma once

// Standard
#include <map>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

// boost
#include <boost/filesystem.hpp>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

namespace fs = boost::filesystem;
namespace dtp {
using namespace seqan3::literals;

using PathVector = std::vector<fs::path>;
using SubPathsMap = std::map<std::string, PathVector>;

using DNAVector = std::vector<seqan3::dna5>;  // seqan3 also provides dna5_vector
using DNASpan = std::span<seqan3::dna5>;
using QualSpan = std::span<seqan3::phred42>;

using QualVector = std::vector<seqan3::phred42>;

// FASTQ
using FASTQFormat = seqan3::type_list<std::string, DNASpan, QualSpan>;
using FASTQFields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>;
using FASTQRecord = seqan3::record<FASTQFormat, FASTQFields>;

// SAM
using sam_field_types =
    seqan3::type_list<std::string, seqan3::sam_flag, std::optional<int32_t>, std::optional<int32_t>,
                      uint8_t, std::vector<seqan3::cigar>, seqan3::dna5_vector,
                      std::vector<seqan3::phred42>, seqan3::sam_tag_dictionary>;

using sam_field_ids =
    seqan3::fields<seqan3::field::id, seqan3::field::flag, seqan3::field::ref_id,
                   seqan3::field::ref_offset, seqan3::field::mapq, seqan3::field::cigar,
                   seqan3::field::seq, seqan3::field::qual, seqan3::field::tags>;

using SamRecord = seqan3::sam_record<sam_field_types, sam_field_ids>;

// Represents a single chimeric read which is split up into its blocks
using SplitRecords = std::vector<SamRecord>;
// Represents a single chimeric read with multiple alternative mapping results
using SplitRecordsVariantGroups = std::vector<SplitRecords>;

// data types used in preprocessing (state transition table)
using Left = std::size_t;                           // left matching block
using Right = std::pair<std::size_t, std::size_t>;  // right matching block

using StateTransitionTable = std::map<std::pair<int, seqan3::dna5>, std::tuple<int, int, int, int>>;
using Tables = std::map<std::pair<std::string, DNAVector>, StateTransitionTable>;
using State = std::tuple<std::string, int, std::pair<int, int>, std::size_t>;

using Bases = std::map<std::pair<std::string, DNAVector>, DNAVector>;
using STTEntry = std::tuple<int, int, int, int>;

enum Strand : char { FORWARD = '+', REVERSE = '-' };

std::optional<int32_t> recordEndPosition(const SamRecord &record);

/**
 * @brief Represents a genomic region.
 */
struct GenomicRegion {
    std::string referenceID;
    int32_t startPosition;
    int32_t endPosition;
    std::optional<Strand> strand;

    /**
     * @brief Constructs a GenomicRegion object.
     * @param referenceID The reference ID of the genomic region.
     * @param startPosition The start position of the genomic region.
     * @param endPosition The end position of the genomic region (exclusive).
     * @param strand The strand of the genomic region (optional).
     */
    GenomicRegion(const std::string &referenceID, int32_t startPosition, int32_t endPosition,
                  std::optional<Strand> strand = std::nullopt);

    /**
     * @brief Creates a GenomicRegion object from a SamRecord.
     * @param record The SamRecord object.
     * @param referenceIDs The deque of reference IDs.
     * @return An optional GenomicRegion object.
     */
    static std::optional<GenomicRegion> fromSamRecord(const SamRecord &record,
                                                      const std::deque<std::string> &referenceIDs);
};

// FeaturesFields for GFF3/GTF
struct Feature {
    std::string referenceID;
    std::string type;
    int32_t startPosition;
    int32_t endPosition;
    Strand strand;
    std::string id;
    std::optional<std::string> groupID;
};

using FeatureMap = std::unordered_map<std::string, std::vector<Feature>>;
}  // namespace dtp

dtp::Strand operator!(dtp::Strand strand);
