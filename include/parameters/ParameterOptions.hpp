#pragma once

// Boost
#include <boost/program_options/options_description.hpp>

// Classes
#include "Constants.hpp"
#include "Orientation.hpp"

namespace po = boost::program_options;

class ParameterOptions {
   public:
    ParameterOptions() = delete;
    ~ParameterOptions() = delete;

    static po::options_description getSubcallOptions();
    static po::options_description getGeneralOptions();
    static po::options_description getPreprocessOptions();
    static po::options_description getAlignOptions();
    static po::options_description getDetectOptions();
    static po::options_description getAnalyzeOptions();
    static po::options_description getOtherOptions();
};
