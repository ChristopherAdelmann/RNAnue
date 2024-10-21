#include "SplitRecordsComplementarityEvaluator.hpp"

auto SplitRecordsComplementarityEvaluator::evaluate(
    const SplitRecords &splitRecords,
    const SplitRecordsEvaluationParameters::BaseParameters &parameters)
    -> std::optional<SplitRecordsComplementarityEvaluator::Result> {
    const auto &sequence1 = splitRecords[0].sequence();
    const auto &sequence2 = splitRecords[1].sequence();

    const auto seq1Reverse = sequence1 | std::views::reverse;

    const auto results =
        CoOptimalPairwiseAligner{SplitRecordsComplementarityEvaluator::complementaryScoringScheme()}
            .getLocalAlignments(std::tie(seq1Reverse, sequence2));

    return SplitRecordsComplementarityEvaluator::getOptimalAlignment(
        results, parameters.minComplementarity, parameters.minComplementarityFraction);
}

constexpr auto SplitRecordsComplementarityEvaluator::complementaryScoringScheme()
    -> seqan3::nucleotide_scoring_scheme<int8_t> {
    using namespace seqan3::literals;

    seqan3::nucleotide_scoring_scheme scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}};

    scheme.score('A'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'A'_dna5) = 1;
    scheme.score('G'_dna5, 'C'_dna5) = 1;
    scheme.score('C'_dna5, 'G'_dna5) = 1;
    scheme.score('G'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'G'_dna5) = 1;

    scheme.score('T'_dna5, 'T'_dna5) = -1;
    scheme.score('A'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'G'_dna5) = -1;

    scheme.score('A'_dna5, 'G'_dna5) = -1;
    scheme.score('A'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'A'_dna5) = -1;

    scheme.score('T'_dna5, 'C'_dna5) = -1;
    scheme.score('C'_dna5, 'T'_dna5) = -1;

    return scheme;
}

auto SplitRecordsComplementarityEvaluator::getOptimalAlignment(
    const std::vector<CoOptimalPairwiseAligner::Result> &alignResults,
    const double minComplementarity, const double minComplementarityFraction)
    -> std::optional<SplitRecordsComplementarityEvaluator::Result> {
    if (alignResults.empty()) {
        return std::nullopt;
    }

    std::optional<CoOptimalPairwiseAligner::Result> optimalAlignment = std::nullopt;

    for (const auto &alignResult : alignResults) {
        if (alignResult.complementarity < minComplementarity ||
            alignResult.fraction < minComplementarityFraction) {
            continue;
        }

        if (!optimalAlignment.has_value() ||
            alignResult.complementarity > optimalAlignment.value().complementarity ||
            (alignResult.complementarity == optimalAlignment.value().complementarity &&
             alignResult.fraction > optimalAlignment.value().fraction)) {
            optimalAlignment = alignResult;
        }
    }

    return optimalAlignment;
}
