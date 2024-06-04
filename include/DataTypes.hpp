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
using SplitRecords = std::vector<SamRecord>;

// data types used in preprocessing (state transition table)
using Left = std::size_t;                           // left matching block
using Right = std::pair<std::size_t, std::size_t>;  // right matching block

using StateTransitionTable = std::map<std::pair<int, seqan3::dna5>, std::tuple<int, int, int, int>>;
using Tables = std::map<std::pair<std::string, DNAVector>, StateTransitionTable>;
using State = std::tuple<std::string, int, std::pair<int, int>, std::size_t>;

using Bases = std::map<std::pair<std::string, DNAVector>, DNAVector>;
using STTEntry = std::tuple<int, int, int, int>;
struct GenomicRegion {
    std::string seqid;
    int start;
    int end;
    GenomicRegion(const std::string &seqid, int start, int end)
        : seqid(seqid), start(start), end(end) {}
};

// FeaturesFields for GFF3/GTF
struct Feature {
    std::string seqid;
    std::string type;
    size_t start;
    size_t end;
    char strand;
    std::string id;
};

using FeatureMap = std::unordered_map<std::string, std::vector<Feature>>;

struct FeatureFields {
    std::string seqid;
    std::string source;
    std::string type;
    int start;
    int end;
    std::string score;
    char strand;
    char phase;
    std::string attributes;
    FeatureFields()
        : seqid(""),
          source(""),
          type(""),
          start(-1),
          end(-1),
          score(""),
          strand(' '),
          phase(' '),
          attributes("") {}
};
}  // namespace dtp
