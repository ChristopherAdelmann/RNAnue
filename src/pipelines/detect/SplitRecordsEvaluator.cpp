#include "SplitRecordsEvaluator.hpp"

SplitRecordsEvaluator::SplitRecordsEvaluator(
    const std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                       SplitRecordsEvaluationParameters::SplicingParameters>
        parameters)
    : parameters(parameters) {}

SplitRecordsEvaluator::Result SplitRecordsEvaluator::evaluate(
    SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs) const {
    if (splitRecords.size() != 2) {
        Logger::log(LogLevel::DEBUG, "Currently only two split records are supported!");
        return SplitRecordsEvaluator::FilterReason::NO_SPLIT_READ;
    }

    if (!splitRecords[0].reference_position().has_value() ||
        !splitRecords[1].reference_position().has_value()) [[unlikely]] {
        Logger::log(LogLevel::WARNING, "Could not determine reference position of split records: ",
                    splitRecords[0].id(), " or ", splitRecords[1].id());
        return SplitRecordsEvaluator::FilterReason::UNMAPPED;
    }

    if (std::holds_alternative<SplitRecordsEvaluationParameters::BaseParameters>(parameters)) {
        return evaluateBase(splitRecords,
                            std::get<SplitRecordsEvaluationParameters::BaseParameters>(parameters));
    } else {
        return evaluateSplicing(
            splitRecords, referenceIDs,
            std::get<SplitRecordsEvaluationParameters::SplicingParameters>(parameters));
    }
}

SplitRecordsEvaluator::Result SplitRecordsEvaluator::evaluateBase(
    SplitRecords &splitRecords,
    const SplitRecordsEvaluationParameters::BaseParameters &parameters) const {
    const auto complementarityResult =
        SplitRecordsComplementarityEvaluator::evaluate(splitRecords, parameters);

    if (!complementarityResult.has_value()) {
        return SplitRecordsEvaluator::FilterReason::COMPLEMENTARITY;
    }

    const auto hybridizationResult =
        SplitRecordsHybridizationEvaluator::evaluate(splitRecords, parameters);

    if (!hybridizationResult.has_value()) {
        return SplitRecordsEvaluator::FilterReason::HYBRIDIZATION;
    }

    addTagsToRecords(splitRecords, complementarityResult.value(), hybridizationResult.value());

    return SplitRecordsEvaluator::EvaluatedSplitRecords{splitRecords, complementarityResult.value(),
                                                        hybridizationResult.value()};
}

SplitRecordsEvaluator::Result SplitRecordsEvaluator::evaluateSplicing(
    SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) const {
    const auto isSplicing =
        SplitRecordsSplicingEvaluator::isSplicedSplitRecord(splitRecords, referenceIDs, parameters);

    if (isSplicing) {
        return SplitRecordsEvaluator::FilterReason::SPLICING;
    }

    return evaluateBase(splitRecords, parameters.baseParameters);
}

void SplitRecordsEvaluator::addTagsToRecords(
    SplitRecords &splitRecords, const SplitRecordsComplementarityEvaluator::Result &complementarity,
    const SplitRecordsHybridizationEvaluator::Result &hybridization) const {
    for (auto &record : splitRecords) {
        // Complementarity tags
        const int length =
            complementarity.endPositions.first - complementarity.beginPositions.first;
        record.tags()["XL"_tag] = length;
        record.tags()["XC"_tag] = static_cast<float>(complementarity.complementarity);
        record.tags()["XR"_tag] = static_cast<float>(complementarity.fraction);
        record.tags()["XS"_tag] = complementarity.score;
        // rec1.tags()["XA"_tag] = res.a;
        // rec2.tags()["XA"_tag] = res.b;
        // rec1.tags()["XM"_tag] = res.matches;
        // rec2.tags()["XM"_tag] = res.matches;

        // Hybridization tags
        record.tags()["XE"_tag] = static_cast<float>(hybridization.energy);
        if (hybridization.crosslinkingResult.has_value()) {
            record.tags()["XN"_tag] =
                static_cast<float>(hybridization.crosslinkingResult.value().normCrosslinkingScore);
            record.tags()["XP"_tag] =
                hybridization.crosslinkingResult.value().preferredCrosslinkingScore;
            record.tags()["XQ"_tag] =
                hybridization.crosslinkingResult.value().nonPreferredCrosslinkingScore;
            record.tags()["XW"_tag] =
                hybridization.crosslinkingResult.value().wobbleCrosslinkingScore;
        }
    }
}

// TODO Implement overall score calculation instead of comparing individual scores
bool SplitRecordsEvaluator::EvaluatedSplitRecords::operator<(
    const EvaluatedSplitRecords &other) const {
    // Calculate if fraction difference is more or equal than 0.05
    if (std::abs(complementarityResult.fraction - other.complementarityResult.fraction) >= 0.05) {
        return complementarityResult.fraction < other.complementarityResult.fraction;
    }

    // Calculate if complementarity difference is more or equal than 0.05
    if (std::abs(complementarityResult.complementarity -
                 other.complementarityResult.complementarity) >= 0.05) {
        return complementarityResult.complementarity < other.complementarityResult.complementarity;
    }

    // Calculate if energy difference is more or equal than 0.05
    if (std::abs(hybridizationResult.energy - other.hybridizationResult.energy) >= 0.05) {
        return hybridizationResult.energy < other.hybridizationResult.energy;
    }

    // If both have crosslinking results prefer the one with higher crosslinking score
    if (hybridizationResult.crosslinkingResult.has_value() &&
        other.hybridizationResult.crosslinkingResult.has_value()) {
        return hybridizationResult.crosslinkingResult.value().normCrosslinkingScore <
               other.hybridizationResult.crosslinkingResult.value().normCrosslinkingScore;
    }

    return false;
}

bool SplitRecordsEvaluator::EvaluatedSplitRecords::operator>(
    const EvaluatedSplitRecords &other) const {
    return !(*this < other);
}

std::ostream &operator<<(std::ostream &os, const SplitRecordsEvaluator::FilterReason &reason) {
    switch (reason) {
        case SplitRecordsEvaluator::FilterReason::NO_SPLIT_READ:
            os << "No split read";
            break;
        case SplitRecordsEvaluator::FilterReason::UNMAPPED:
            os << "Unmapped";
            break;
        case SplitRecordsEvaluator::FilterReason::SPLICING:
            os << "Splicing";
            break;
        case SplitRecordsEvaluator::FilterReason::COMPLEMENTARITY:
            os << "Complementarity";
            break;
        case SplitRecordsEvaluator::FilterReason::HYBRIDIZATION:
            os << "Hybridization";
            break;
    }

    return os;
}
