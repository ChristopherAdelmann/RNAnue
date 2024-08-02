#pragma once

// Standard
#include <climits>
#include <cstddef>
#include <filesystem>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Classes
#include "GeneralParameters.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;

class AlignParameters : public GeneralParameters {
   public:
    std::filesystem::path referenceGenome;
    size_t accuracy;
    size_t minimumFragmentScore;
    size_t minimumFragmentLength;
    size_t minimumSpliceCoverage;

    AlignParameters(const po::variables_map& params)
        : GeneralParameters(params),
          referenceGenome(ParameterValidator::validateFilePath(params, "dbref")),
          accuracy(ParameterValidator::validateArithmetic(params, "accuracy", 0, 100)),
          minimumFragmentScore(
              ParameterValidator::validateArithmetic(params, "minfragsco", 0, INT_MAX)),
          minimumFragmentLength(
              ParameterValidator::validateArithmetic(params, "minfraglen", 0, INT_MAX)),
          minimumSpliceCoverage(
              ParameterValidator::validateArithmetic(params, "minsplicecov", 0, 100)) {};
};
