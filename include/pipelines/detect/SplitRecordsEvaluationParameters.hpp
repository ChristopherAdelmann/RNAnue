#pragma once

// Standard
#include <cstdint>
#include <variant>

// Classes
#include "FeatureAnnotator.hpp"
#include "Orientation.hpp"

namespace SplitRecordsEvaluationParameters {
struct BaseParameters {
    double minComplementarity;
    double minComplementarityFraction;
    double mfeThreshold;
};

struct SplicingParameters {
    BaseParameters baseParameters;
    Annotation::Orientation orientation;
    int32_t splicingTolerance;
    const Annotation::FeatureAnnotator& featureAnnotator;
};

using ParameterVariant = std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                                      SplitRecordsEvaluationParameters::SplicingParameters>;
}  // namespace SplitRecordsEvaluationParameters
