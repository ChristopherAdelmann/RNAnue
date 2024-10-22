#pragma once

// Standard
#include <cfloat>
#include <climits>
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "GeneralParameters.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;

namespace pipelines::detect {

struct DetectParameters : public GeneralParameters {
   public:
    size_t minimumFragmentLength;
    size_t minimumMapQuality;
    double minimumComplementarity;
    double minimumSiteLengthRatio;
    double maxHybridizationEnergy;
    bool excludeSoftClipping;
    bool removeSplicingEvents;
    int splicingTolerance;

    DetectParameters(const po::variables_map& params)
        : GeneralParameters(params),
          minimumFragmentLength(
              ParameterValidator::validateArithmetic<size_t>(params, "minfraglen", 0, SIZE_MAX)),
          minimumMapQuality(
              ParameterValidator::validateArithmetic<size_t>(params, "mapqmin", 0, SIZE_MAX)),
          minimumComplementarity(
              ParameterValidator::validateArithmetic(params, "cmplmin", 0.0, 1.0)),
          minimumSiteLengthRatio(
              ParameterValidator::validateArithmetic(params, "sitelenratio", 0.0, 1.0)),
          maxHybridizationEnergy(
              ParameterValidator::validateArithmetic(params, "nrgmax", DBL_MIN, DBL_MAX)),
          excludeSoftClipping(validateExcludeSoftClipping(params)),
          removeSplicingEvents(validateRemoveSplicingEvents(params)),
          splicingTolerance(ParameterValidator::validateArithmetic(params, "splicingtolerance",
                                                                   INT_MIN, INT_MAX)) {};

   private:
    static auto validateExcludeSoftClipping(const po::variables_map& params) -> bool {
        return params["exclclipping"].as<bool>();
    }

    static auto validateRemoveSplicingEvents(const po::variables_map& params) -> bool {
        return params["splicing"].as<bool>();
    }
};

}  // namespace pipelines::detect
