#pragma once

// Standard
#include <climits>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <istream>
#include <set>
#include <string>
#include <unordered_set>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Classes
#include "Logger.hpp"
#include "Orientation.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;
class GeneralParameters {
   public:
    std::filesystem::path treatmentsDir;
    std::optional<std::filesystem::path> controlDir;
    std::filesystem::path outputDir;

    std::filesystem::path featuresInPath;
    std::unordered_set<std::string> featureTypes;
    Annotation::Orientation featureOrientation;

    LogLevel logLevel;

    size_t threadCount;
    size_t chunkSize;

    GeneralParameters(const po::variables_map& params)
        : treatmentsDir(ParameterValidator::validateDirectory(params, "trtms")),
          controlDir(validateControlDir(params)),
          outputDir(ParameterValidator::validateDirectory(params, "outdir")),
          featuresInPath(ParameterValidator::validateFilePath(params, "features")),
          featureTypes(validateFeatureTypes(params)),
          featureOrientation(validateFeatureOrientation(params)),
          logLevel(validateLogLevel(params)),
          threadCount(ParameterValidator::validateArithmetic(params, "threads", 1, INT_MAX)),
          chunkSize(ParameterValidator::validateArithmetic(params, "chunksize", 1, INT_MAX)) {};

   private:
    static std::optional<std::filesystem::path> validateControlDir(
        const po::variables_map& params) {
        if (params.count("ctrls") && !params["ctrls"].as<std::string>().empty()) {
            return ParameterValidator::validateDirectory(params, "ctrls");
        } else {
            return std::nullopt;
        }
    }

    static std::unordered_set<std::string> validateFeatureTypes(const po::variables_map& params) {
        const auto featureTypesString = params["featuretypes"].as<std::string>();

        std::unordered_set<std::string> uniqueIncludedFeatures;

        std::stringstream ss(featureTypesString);
        std::string str;
        while (getline(ss, str, ',')) {
            str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
            uniqueIncludedFeatures.insert(str);
        }

        return uniqueIncludedFeatures;
    }

    static Annotation::Orientation validateFeatureOrientation(const po::variables_map& params) {
        return params["orientation"].as<Annotation::Orientation>();
    }

    static LogLevel validateLogLevel(const po::variables_map& params) {
        const std::string logLevelStr = params["loglevel"].as<std::string>();

        if (logLevelStr == "debug" || logLevelStr == "DEBUG") {
            return LogLevel::DEBUG;
        } else if (logLevelStr == "info" || logLevelStr == "INFO") {
            return LogLevel::INFO;
        } else if (logLevelStr == "warning" || logLevelStr == "WARNING") {
            return LogLevel::WARNING;
        } else if (logLevelStr == "error" || logLevelStr == "ERROR") {
            return LogLevel::ERROR;
        } else {
            Logger::log(LogLevel::ERROR, "Invalid log level specified.");
            return LogLevel::INFO;
        }
    }
};
