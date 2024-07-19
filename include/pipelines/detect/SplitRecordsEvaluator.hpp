#pragma once

// Standard
#include <deque>
#include <iostream>
#include <optional>
#include <variant>

// Classes
#include "CooptimalPairwiseAligner.hpp"
#include "DataTypes.hpp"
#include "FeatureAnnotator.hpp"
#include "SplitRecordsComplementarityEvaluator.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsHybridizationEvaluator.hpp"
#include "SplitRecordsSplicingEvaluator.hpp"

using namespace seqan3::literals;

class SplitRecordsEvaluator {
   public:
    struct EvaluatedSplitRecords {
        dtp::SplitRecords splitRecords;
        SplitRecordsComplementarityEvaluator::Result complementarityResult;
        SplitRecordsHybridizationEvaluator::Result hybridizationResult;

        bool operator<(const EvaluatedSplitRecords &other) const;
        bool operator>(const EvaluatedSplitRecords &other) const;
    };

    enum class FilterReason { NO_SPLIT_READ, UNMAPPED, SPLICING, COMPLEMENTARITY, HYBRIDIZATION };

    using Result = std::variant<EvaluatedSplitRecords, FilterReason>;

    explicit SplitRecordsEvaluator(
        const std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                           SplitRecordsEvaluationParameters::SplicingParameters>
            parameters);
    ~SplitRecordsEvaluator() = default;

    Result evaluate(dtp::SplitRecords &splitRecords,
                    const std::deque<std::string> &referenceIDs) const;

   private:
    std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                 SplitRecordsEvaluationParameters::SplicingParameters>
        parameters;

    Result evaluateBase(dtp::SplitRecords &splitRecords,
                        const SplitRecordsEvaluationParameters::BaseParameters &parameters) const;

    Result evaluateSplicing(
        dtp::SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters) const;

    void addTagsToRecords(dtp::SplitRecords &splitRecords,
                          const SplitRecordsComplementarityEvaluator::Result &complementarity,
                          const SplitRecordsHybridizationEvaluator::Result &hybridization) const;
};

std::ostream &operator<<(std::ostream &os, const SplitRecordsEvaluator::FilterReason &reason);