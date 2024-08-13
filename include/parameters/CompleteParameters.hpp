#pragma once

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Classes
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "DetectParameters.hpp"
#include "PreprocessParameters.hpp"

namespace po = boost::program_options;

struct CompleteParameters {
    PreprocessParameters preprocessParameters;
    AlignParameters alignParameters;
    DetectParameters detectParameters;
    AnalyzeParameters analyzeParameters;

    CompleteParameters(const po::variables_map &params)
        : preprocessParameters(params),
          alignParameters(params),
          detectParameters(params),
          analyzeParameters(params) {};
};
