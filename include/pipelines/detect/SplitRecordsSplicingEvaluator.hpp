#pragma once

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

using namespace dataTypes;

class SplitRecordsSplicingEvaluator {
   public:
    SplitRecordsSplicingEvaluator() = delete;
    ~SplitRecordsSplicingEvaluator() = delete;
    SplitRecordsSplicingEvaluator(const SplitRecordsSplicingEvaluator &) = delete;
    auto operator=(const SplitRecordsSplicingEvaluator &) -> SplitRecordsSplicingEvaluator & =
                                                                 delete;
    SplitRecordsSplicingEvaluator(SplitRecordsSplicingEvaluator &&) = delete;
    auto operator=(SplitRecordsSplicingEvaluator &&) -> SplitRecordsSplicingEvaluator & = delete;

    static auto isSplicedSplitRecord(
        const SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> bool;

   private:
    static auto getGroupedFeatures(
        const SamRecord &record1, const SamRecord &record2,
        const std::deque<std::string> &referenceIDs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters)
        -> std::optional<std::pair<dataTypes::GenomicFeature, dataTypes::GenomicFeature>>;

    static auto getSpliceJunctionBoundingRegion(
        const SamRecord &record1, const SamRecord &record2,
        const dataTypes::GenomicFeature &featureRecord1,
        const dataTypes::GenomicFeature &featureRecord2,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters)
        -> std::optional<dataTypes::GenomicRegion>;
};
