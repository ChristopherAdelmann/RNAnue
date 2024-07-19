#pragma once

// Classes
#include "DataTypes.hpp"
#include "FeatureAnnotator.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "Utility.hpp"

class SplitRecordsSplicingEvaluator {
   public:
    SplitRecordsSplicingEvaluator() = delete;
    ~SplitRecordsSplicingEvaluator() = delete;
    SplitRecordsSplicingEvaluator(const SplitRecordsSplicingEvaluator &) = delete;
    SplitRecordsSplicingEvaluator &operator=(const SplitRecordsSplicingEvaluator &) = delete;
    SplitRecordsSplicingEvaluator(SplitRecordsSplicingEvaluator &&) = delete;
    SplitRecordsSplicingEvaluator &operator=(SplitRecordsSplicingEvaluator &&) = delete;

    static bool isSplicedSplitRecord(
        const dtp::SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters);

   private:
    static std::optional<std::pair<dtp::Feature, dtp::Feature>> getGroupedFeatures(
        const dtp::SamRecord &record1, const dtp::SamRecord &record2,
        const std::deque<std::string> &referenceIDs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters);

    static std::optional<dtp::GenomicRegion> getSpliceJunctionBoundingRegion(
        const dtp::SamRecord &record1, const dtp::SamRecord &record2,
        const dtp::Feature &featureRecord1, const dtp::Feature &featureRecord2,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters);
};