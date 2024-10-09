#pragma once

// Standard
#include <map>
#include <optional>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

// ViennaRNA
extern "C" {
#include <ViennaRNA/cofold.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
}

// Classes
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

using namespace seqan3::literals;
using namespace dataTypes;

class SplitRecordsHybridizationEvaluator {
   public:
    struct Result {
        struct CrosslinkingResult {
            double normCrosslinkingScore;
            int preferredCrosslinkingScore;
            int nonPreferredCrosslinkingScore;
            int wobbleCrosslinkingScore;
        };

        double energy;
        std::optional<CrosslinkingResult> crosslinkingResult;
    };

    SplitRecordsHybridizationEvaluator() = delete;
    ~SplitRecordsHybridizationEvaluator() = delete;
    SplitRecordsHybridizationEvaluator(const SplitRecordsHybridizationEvaluator &) = delete;
    SplitRecordsHybridizationEvaluator &operator=(const SplitRecordsHybridizationEvaluator &) =
        delete;
    SplitRecordsHybridizationEvaluator(SplitRecordsHybridizationEvaluator &&) = delete;
    SplitRecordsHybridizationEvaluator &operator=(SplitRecordsHybridizationEvaluator &&) = delete;

    static std::optional<Result> evaluate(
        const SplitRecords &splitRecords,
        const SplitRecordsEvaluationParameters::BaseParameters &parameters);

   private:
    using NucleotidePairPositions = std::pair<size_t, size_t>;
    using NucleotidePositionsWindow = std::pair<NucleotidePairPositions, NucleotidePairPositions>;
    using NucleotideWindowPair = std::pair<seqan3::dna5_vector, seqan3::dna5_vector>;
    static const std::map<NucleotideWindowPair, size_t> crosslinkingScoringScheme;

    struct InteractionWindow {
        seqan3::dna5_vector forwardWindowNucleotides;
        seqan3::dna5_vector reverseWindowNucleotides;
        std::pair<size_t, size_t> forwardWindowPositions;
        std::pair<size_t, size_t> reverseWindowPositions;
        bool isInterFragment;
    };

    static std::optional<Result::CrosslinkingResult> findCrosslinkingSites(
        std::span<const seqan3::dna5> sequence1, std::span<const seqan3::dna5> sequence2,
        std::vector<seqan3::dot_bracket3> &dotbracket);

    static std::optional<InteractionWindow> getContinuosNucleotideWindows(
        std::span<const seqan3::dna5> sequence1, std::span<const seqan3::dna5> sequence2,
        NucleotidePositionsWindow positionsPair);
};
