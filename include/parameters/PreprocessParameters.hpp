#pragma once

// Standard
#include <climits>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <string>
#include <variant>
#include <vector>

// Boost
#include <boost/program_options/variables_map.hpp>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/range/to.hpp>

// Classes
#include "Constants.hpp"
#include "GeneralParameters.hpp"
#include "Logger.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;

class PreprocessParameters : public GeneralParameters {
   public:
    bool trimPolyG;
    bool preprocessEnabled;

    using Adapter = std::variant<std::filesystem::path, seqan3::dna5_vector>;
    std::optional<Adapter> adapter5Forward;
    std::optional<Adapter> adapter3Forward;
    std::optional<Adapter> adapter5Reverse;
    std::optional<Adapter> adapter3Reverse;

    double maxMissMatchFractionTrimming;
    size_t minOverlapTrimming;
    size_t minQualityThreshold;
    size_t minLengthThreshold;
    size_t minMeanWindowQuality;
    size_t windowTrimmingSize;
    size_t minOverlapMerging;
    double maxMissMatchFractionMerging;

    PreprocessParameters(const po::variables_map& params)
        : GeneralParameters(params),
          trimPolyG(validateTrimPolyG(params)),
          preprocessEnabled(validatePreprocessEnabled(params)),
          adapter5Forward(validateAdapter(params, "adapter5Forward")),
          adapter3Forward(validateAdapter(params, "adapter3Forward")),
          adapter5Reverse(validateAdapter(params, "adapter5Reverse")),
          adapter3Reverse(validateAdapter(params, "adapter3Reverse")),
          maxMissMatchFractionTrimming(
              ParameterValidator::validateArithmetic(params, "mtrim", 0.0, 1.0)),
          minOverlapTrimming(
              ParameterValidator::validateArithmetic(params, "minovltrim", size_t{0}, SIZE_T_MAX)),
          minQualityThreshold(
              ParameterValidator::validateArithmetic(params, "minqual", size_t{0}, SIZE_T_MAX)),
          minLengthThreshold(
              ParameterValidator::validateArithmetic(params, "minlen", size_t{0}, SIZE_T_MAX)),
          minMeanWindowQuality(
              ParameterValidator::validateArithmetic(params, "wqual", size_t{0}, SIZE_T_MAX)),
          windowTrimmingSize(
              ParameterValidator::validateArithmetic(params, "wtrim", size_t{0}, SIZE_T_MAX)),
          minOverlapMerging(
              ParameterValidator::validateArithmetic(params, "minovl", size_t{0}, SIZE_T_MAX)),
          maxMissMatchFractionMerging(
              ParameterValidator::validateArithmetic(params, "mmerge", 0.0, 1.0)) {}

   private:
    static bool validateTrimPolyG(const po::variables_map& params) {
        return params["trimpolyg"].as<bool>();
    }

    static bool validatePreprocessEnabled(const po::variables_map& params) {
        return params[constants::pipelines::PREPROCESS.c_str()].as<bool>();
    }

    static std::optional<Adapter> validateAdapter(const po::variables_map& params,
                                                  const std::string& paramName) {
        if (params.count(paramName)) {
            const std::string adapterStr = params[paramName].as<std::string>();
            if (std::filesystem::exists(adapterStr)) {
                return std::filesystem::path(adapterStr);
            } else {
                return adapterStr | seqan3::views::char_to<seqan3::dna5> |
                       seqan3::ranges::to<std::vector>();
            }
        }
        return std::nullopt;
    }
};