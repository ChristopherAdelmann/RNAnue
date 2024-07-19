#pragma once

#include "Annotation/FeatureAnnotator.hpp"

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
    const Annotation::FeatureAnnotator &featureAnnotator;
};
}