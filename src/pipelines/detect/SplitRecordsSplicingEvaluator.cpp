#include "SplitRecordsSplicingEvaluator.hpp"

bool SplitRecordsSplicingEvaluator::isSplicedSplitRecord(
    const SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) {
    if (splitRecords[0].reference_id() != splitRecords[1].reference_id()) {
        return false;
    }

    const auto &record1 =
        splitRecords[0].reference_position().value() < splitRecords[1].reference_position().value()
            ? splitRecords[0]
            : splitRecords[1];
    const auto &record2 =
        splitRecords[0].reference_position().value() < splitRecords[1].reference_position().value()
            ? splitRecords[1]
            : splitRecords[0];

    const auto features = getGroupedFeatures(record1, record2, referenceIDs, parameters);

    if (!features.has_value()) {
        return false;
    }

    const auto &featureRecord1 = features.value().first;
    const auto &featureRecord2 = features.value().second;

    const auto spliceJunctionBoundingRegion = getSpliceJunctionBoundingRegion(
        record1, record2, featureRecord1, featureRecord2, parameters);

    if (!spliceJunctionBoundingRegion.has_value()) {
        return false;
    }

    auto it = parameters.featureAnnotator.overlappingFeatureIterator(
        spliceJunctionBoundingRegion.value(), parameters.orientation);

    const bool hasInBetweenExon =
        std::any_of(it.begin(), it.end(), [&featureRecord1](const auto &feature) {
            return feature.groupID.has_value() && feature.groupID.value() == featureRecord1.groupID;
        });

    if (hasInBetweenExon) {
        return false;
    }

    return true;
}

std::optional<std::pair<dtp::Feature, dtp::Feature>>
SplitRecordsSplicingEvaluator::getGroupedFeatures(
    const SamRecord &record1, const SamRecord &record2, const std::deque<std::string> &referenceIDs,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) {
    const auto featureRecord1 = parameters.featureAnnotator.getBestOverlappingFeature(
        record1, referenceIDs, parameters.orientation);

    if (!featureRecord1.has_value() || !featureRecord1.value().groupID.has_value()) {
        return std::nullopt;
    }

    const auto featureRecord2 = parameters.featureAnnotator.getBestOverlappingFeature(
        record2, referenceIDs, parameters.orientation);

    if (!featureRecord2.has_value() || !featureRecord2.value().groupID.has_value()) {
        return std::nullopt;
    }

    if (featureRecord1.value().groupID != featureRecord2.value().groupID) {
        return std::nullopt;
    }

    return std::make_pair(featureRecord1.value(), featureRecord2.value());
}

std::optional<dtp::GenomicRegion> SplitRecordsSplicingEvaluator::getSpliceJunctionBoundingRegion(
    const SamRecord &record1, const SamRecord &record2, const dtp::Feature &featureRecord1,
    const dtp::Feature &featureRecord2,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) {
    const auto record1EndPosition = dtp::recordEndPosition(record1).value();

    bool record1AtSpliceJunctionStart = helper::isContained(
        record1EndPosition, featureRecord1.endPosition - 1, parameters.splicingTolerance);

    if (!record1AtSpliceJunctionStart) {
        return std::nullopt;
    }

    const auto record2StartPosition = record2.reference_position().value();
    bool record2AtSpliceJunctionEnd = helper::isContained(
        record2StartPosition, featureRecord2.startPosition - 1, parameters.splicingTolerance);

    if (!record2AtSpliceJunctionEnd) {
        return std::nullopt;
    }

    return dtp::GenomicRegion{featureRecord1.referenceID, record1EndPosition + 1,
                              record2StartPosition};
}
