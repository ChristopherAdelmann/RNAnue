#pragma once

// Standard
#include <ranges>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/utility/range/to.hpp>

// varna
extern "C" {
#include <ViennaRNA/cofold.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
}

// Classes
#include "CooptimalPairwiseAligner.hpp"
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"
#include "Logger.hpp"

using namespace dtp;
using namespace seqan3::literals;

using NucleotidePairPositions = std::pair<size_t, size_t>;
using NucleotidePositionsWindow = std::pair<NucleotidePairPositions, NucleotidePairPositions>;
using NucleotideWindowPair = std::pair<seqan3::dna5_vector, seqan3::dna5_vector>;

static const std::map<NucleotideWindowPair, size_t> crosslinkingScoringScheme = {
    {{"TA"_dna5, "AT"_dna5}, 3},  // Preferred pyrimidine cross-linking
    {{"TG"_dna5, "GT"_dna5}, 3},
    {{"TC"_dna5, "CT"_dna5}, 3},
    {{"AT"_dna5, "TA"_dna5}, 3},
    {{"GT"_dna5, "TG"_dna5}, 3},
    {{"CT"_dna5, "TC"_dna5}, 3},
    {{"CA"_dna5, "AC"_dna5}, 2},  // Non-preferred pyrimidine cross-linking
    {{"CG"_dna5, "GC"_dna5}, 2},
    {{"AC"_dna5, "CA"_dna5}, 2},
    {{"GC"_dna5, "CG"_dna5}, 2},
    {{"TA"_dna5, "GT"_dna5}, 1},  // Wobble base pairs cross-linking
    {{"TG"_dna5, "GT"_dna5}, 1},
    {{"GT"_dna5, "TA"_dna5}, 1},
    {{"GT"_dna5, "TG"_dna5}, 1}};

struct InteractionWindow {
    seqan3::dna5_vector forwardWindowNucleotides;
    seqan3::dna5_vector reverseWindowNucleotides;
    std::pair<size_t, size_t> forwardWindowPositions;
    std::pair<size_t, size_t> reverseWindowPositions;
    bool isInterFragment;
};

struct CrosslinkingResult {
    double normCrosslinkingScore;
    int preferredCrosslinkingScore;
    int nonPreferredCrosslinkingScore;
    int wobbleCrosslinkingScore;
};

struct HybridizationResult {
    double energy;
    std::optional<CrosslinkingResult> crosslinkingResult;
};

class EvaluatedSplitRecords {
   public:
       static std::optional<EvaluatedSplitRecords> calculateEvaluatedSplitRecords(
        SplitRecords &splitRecords, const double minComplementarity,
        const double minComplementarityFraction, const double mfeThreshold);

    const SplitRecords splitRecords;
    const CoOptimalPairwiseAligner::AlignmentResult complementarityResult;
    const HybridizationResult hybridizationResult;

    bool operator<(const EvaluatedSplitRecords &other) const;
    bool operator>(const EvaluatedSplitRecords &other) const;

   private:
    EvaluatedSplitRecords(SplitRecords splitRecords,
                          CoOptimalPairwiseAligner::AlignmentResult complementarity,
                          HybridizationResult hybridization)
        : splitRecords(splitRecords),
          complementarityResult(complementarity),
          hybridizationResult(hybridization) {}

    // Add the complementarity and hybridization results to the split records tags
    static void addTagsToRecords(SplitRecords &splitRecords,
                                 const CoOptimalPairwiseAligner::AlignmentResult &complementarity,
                                 const HybridizationResult &hybridization);

    // Calculate the complementarity between two sequences
    static std::optional<CoOptimalPairwiseAligner::AlignmentResult> calculateComplementarity(
        const seqan3::dna5_vector &seq1, const seqan3::dna5_vector &seq2,
        const double minComplementarity, const double minComplementarityFraction);

    static constexpr seqan3::nucleotide_scoring_scheme<int8_t> complementaryScoringScheme();

    static std::optional<CoOptimalPairwiseAligner::AlignmentResult> getOptimalAlignment(
        const std::vector<CoOptimalPairwiseAligner::AlignmentResult> &alignResults,
        const double minComplementarity, const double minComplementarityFraction);

    // Calculate the hybridization between two sequences
    static std::optional<HybridizationResult> calculateHybridization(
        std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
        const double mfeThreshold);

    static std::optional<CrosslinkingResult> findCrosslinkingSites(
        std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
        std::vector<seqan3::dot_bracket3> &dotbracket);
    static std::optional<InteractionWindow> getContinuosNucleotideWindows(
        std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
        NucleotidePositionsWindow positionsPair);
};