#pragma once

// Boost
#include <boost/program_options/options_description.hpp>

namespace po = boost::program_options;

class ParameterOptions {
   public:
    ParameterOptions() = delete;
    ParameterOptions(const ParameterOptions &) = default;
    ParameterOptions(ParameterOptions &&) = delete;
    auto operator=(const ParameterOptions &) -> ParameterOptions & = default;
    auto operator=(ParameterOptions &&) -> ParameterOptions & = delete;
    ~ParameterOptions() = delete;

    static auto getSubcallOptions() -> po::options_description;
    static auto getGeneralOptions() -> po::options_description;
    static auto getPreprocessOptions() -> po::options_description;
    static auto getAlignOptions() -> po::options_description;
    static auto getDetectOptions() -> po::options_description;
    static auto getAnalyzeOptions() -> po::options_description;
    static auto getOtherOptions() -> po::options_description;
};
