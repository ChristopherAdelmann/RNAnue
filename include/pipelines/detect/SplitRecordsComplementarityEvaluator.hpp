#pragma once

// Standard
#include <optional>
#include <variant>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Classes
#include "CooptimalPairwiseAligner.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

class SplitRecordsComplementarityEvaluator {
   public:
    using Result = CoOptimalPairwiseAligner::Result;

    SplitRecordsComplementarityEvaluator() = delete;
    ~SplitRecordsComplementarityEvaluator() = delete;
    SplitRecordsComplementarityEvaluator(const SplitRecordsComplementarityEvaluator &) = delete;
    SplitRecordsComplementarityEvaluator &operator=(const SplitRecordsComplementarityEvaluator &) =
        delete;
    SplitRecordsComplementarityEvaluator(SplitRecordsComplementarityEvaluator &&) = delete;
    SplitRecordsComplementarityEvaluator &operator=(SplitRecordsComplementarityEvaluator &&) =
        delete;

    static std::optional<SplitRecordsComplementarityEvaluator::Result> evaluate(
        const SplitRecords &splitRecords,
        const SplitRecordsEvaluationParameters::BaseParameters &parameters);

   private:
    static constexpr seqan3::nucleotide_scoring_scheme<int8_t> complementaryScoringScheme();

    static std::optional<SplitRecordsComplementarityEvaluator::Result> getOptimalAlignment(
        const std::vector<CoOptimalPairwiseAligner::Result> &alignResults,
        const double minComplementarity, const double minComplementarityFraction);
};
