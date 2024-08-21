#pragma once

#include <optional>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "seqan3/alignment/configuration/align_config_gap_cost_affine.hpp"
#include "seqan3/alignment/configuration/align_config_method.hpp"
#include "seqan3/alignment/configuration/align_config_output.hpp"
#include "seqan3/alignment/configuration/align_config_scoring_scheme.hpp"
#include "seqan3/alignment/pairwise/align_pairwise.hpp"
#include "seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp"
#include "seqan3/alphabet/quality/phred42.hpp"
#include "seqan3/alphabet/views/complement.hpp"
namespace pipelines {
namespace preprocess {

struct PairedRecordMerger {
    /**
     * @brief Merges two records if they have a sufficient overlap.
     *
     * This function takes two records and checks if they have a sufficient overlap based on the
     * specified minimum overlap value. If the overlap is sufficient, the records are merged into a
     * new record and returned as an optional. Otherwise, std::nullopt is returned.
     *
     * @tparam record_type The type of the records.
     * @param record1 The first record to be merged.
     * @param record2 The second record to be merged.
     * @param minOverlap The minimum required overlap between the records.
     * @return An optional containing the merged record if the overlap is sufficient, otherwise
     * std::nullopt.
     */
    template <typename record_type>
    static std::optional<record_type> mergeRecordPair(const record_type &record1,
                                                      const record_type &record2,
                                                      const size_t minOverlapMerge,
                                                      const double maxMissMatchRateMerge) {
        const seqan3::align_cfg::method_global endGapConfig{
            seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
            seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};

        const seqan3::align_cfg::scoring_scheme scoringSchemeConfig{
            seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};

        const seqan3::align_cfg::gap_cost_affine gapSchemeConfig{
            seqan3::align_cfg::open_score{-2}, seqan3::align_cfg::extension_score{-4}};

        const auto outputConfig =
            seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{} |
            seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

        const auto alignmentConfig =
            endGapConfig | scoringSchemeConfig | gapSchemeConfig | outputConfig;

        const auto &seq1 = record1.sequence();
        const auto &seq2ReverseComplement =
            record2.sequence() | std::views::reverse | seqan3::views::complement;

        std::optional<record_type> mergedRecord{std::nullopt};

        for (auto const &result :
             seqan3::align_pairwise(std::tie(seq1, seq2ReverseComplement), alignmentConfig)) {
            const int overlap = result.sequence1_end_position() - result.sequence1_begin_position();
            if (overlap < int(minOverlapMerge)) continue;

            const int minScore = overlap - (overlap * maxMissMatchRateMerge) * 2;
            if (result.score() >= minScore) {
                mergedRecord = constructMergedRecord(record1, record2, result);
            }
        }

        return mergedRecord;
    }

    /**
     * Constructs a merged record by combining two input records with a specified overlap.
     *
     * @tparam record_type The type of the input records.
     * @param record1 The first input record.
     * @param record2 The second input record.
     * @param overlap The length of the overlap between the two records.
     * @return The merged record.
     */
    template <typename record_type, typename result_type>
    static record_type constructMergedRecord(
        const record_type &record1, const record_type &record2,
        const seqan3::alignment_result<result_type> &alignmentResult) {
        const auto &record1Qualities = record1.base_qualities();
        const auto &record2Qualities = record2.base_qualities();

        seqan3::dna5_vector mergedSequence{};
        std::vector<seqan3::phred42> mergedQualities{};

        mergedSequence.reserve(record1.sequence().size() + record2.sequence().size());
        mergedQualities.reserve(record1Qualities.size() + record2Qualities.size());

        const auto &record2RevCompSequence =
            record2.sequence() | seqan3::views::complement | std::views::reverse;
        const auto &record2RevCompQualities = record2Qualities | std::views::reverse;

        // 5' overhang of read one that is not in overlap region
        mergedSequence.insert(
            mergedSequence.end(), record1.sequence().begin(),
            record1.sequence().begin() + alignmentResult.sequence1_begin_position());
        mergedQualities.insert(
            mergedQualities.end(), record1.base_qualities().begin(),
            record1.base_qualities().begin() + alignmentResult.sequence1_begin_position());

        size_t posRecord1 = alignmentResult.sequence1_begin_position();
        size_t posRecord2 = alignmentResult.sequence2_begin_position();

        const auto &[alignmentSeq1, alignmentSeq2] = alignmentResult.alignment();

        for (const auto &[el1, el2] : seqan3::views::zip(alignmentSeq1, alignmentSeq2)) {
            const seqan3::phred42 qual1 = record1Qualities[posRecord1];
            const seqan3::phred42 qual2 = record2RevCompQualities[posRecord2];

            if (el1 == el2) {
                const seqan3::dna5 base1 = el1.template convert_to<seqan3::dna5>();
                mergedSequence.insert(mergedSequence.end(), base1);
                mergedQualities.insert(mergedQualities.end(), std::max(qual1, qual2));

                posRecord1++;
                posRecord2++;
                continue;
            }

            if (el1 != seqan3::gap{} && el2 != seqan3::gap{}) {
                if (qual1 < qual2) {
                    const seqan3::dna5 base2 = el2.template convert_to<seqan3::dna5>();
                    mergedSequence.insert(mergedSequence.end(), base2);
                    mergedQualities.insert(mergedQualities.end(), qual2);
                } else {
                    const seqan3::dna5 base1 = el1.template convert_to<seqan3::dna5>();
                    mergedSequence.insert(mergedSequence.end(), base1);
                    mergedQualities.insert(mergedQualities.end(), qual1);
                }

                posRecord1++;
                posRecord2++;
                continue;
            }

            if (el1 == seqan3::gap{}) {
                const seqan3::dna5 base2 = el2.template convert_to<seqan3::dna5>();
                mergedSequence.insert(mergedSequence.end(), base2);
                mergedQualities.insert(mergedQualities.end(), qual2);

                posRecord2++;
                continue;
            }

            if (el2 == seqan3::gap{}) {
                const seqan3::dna5 base1 = el1.template convert_to<seqan3::dna5>();
                mergedSequence.insert(mergedSequence.end(), base1);
                mergedQualities.insert(mergedQualities.end(), qual1);

                posRecord1++;
                continue;
            }
        }

        // 5' overhang of read two that is not in overlap region
        mergedSequence.insert(
            mergedSequence.end(),
            record2RevCompSequence.begin() + alignmentResult.sequence2_end_position(),
            record2RevCompSequence.end());
        mergedQualities.insert(
            mergedQualities.end(),
            record2RevCompQualities.begin() + alignmentResult.sequence2_end_position(),
            record2RevCompQualities.end());

        return record_type{std::move(mergedSequence), record1.id(), std::move(mergedQualities)};
    }
};

}  // namespace preprocess
}  // namespace pipelines
