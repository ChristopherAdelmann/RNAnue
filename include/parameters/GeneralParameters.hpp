#pragma once

// Standard
#include <climits>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <istream>
#include <set>
#include <string>
#include <vector>

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
    enum class ReadType { SINGLE_END, PAIRED_END };

    // Parameters
    ReadType readType;

    std::filesystem::path treatmentsDir;
    std::filesystem::path controlDir;
    std::filesystem::path outputDir;

    std::filesystem::path featuresInPath;
    std::vector<std::string> featureTypes;
    Annotation::Orientation featureOrientation;

    LogLevel logLevel;

    size_t threadCount;
    size_t chunkSize;

    GeneralParameters(const po::variables_map& params)
        : readType(validateReadType(params)),
          treatmentsDir(ParameterValidator::validateDirectory(params, "trtms")),
          controlDir(ParameterValidator::validateDirectory(params, "ctrls")),
          outputDir(ParameterValidator::validateDirectory(params, "outdir")),
          featuresInPath(ParameterValidator::validateFilePath(params, "features")),
          featureTypes(validateFeatureTypes(params)),
          featureOrientation(validateFeatureOrientation(params)),
          logLevel(validateLogLevel(params)),
          threadCount(ParameterValidator::validateArithmetic(params, "threads", 1, INT_MAX)),
          chunkSize(ParameterValidator::validateArithmetic(params, "chunksize", 1, INT_MAX)) {};

   private:
    static ReadType validateReadType(const po::variables_map& params) {
        const std::string readType = params["readtype"].as<std::string>();
        const std::set<std::string> availableReadTypes{"SE", "PE"};

        if (!availableReadTypes.contains(readType)) {
            Logger::log(LogLevel::ERROR,
                        "Parameter readtype is none of the expected values [SE, PE].");
        }

        return readType == "SE" ? ReadType::SINGLE_END : ReadType::PAIRED_END;
    }

    static std::vector<std::string> validateFeatureTypes(const po::variables_map& params) {
        return params["featuretypes"].as<std::vector<std::string>>();
    }

    static Annotation::Orientation validateFeatureOrientation(const po::variables_map& params) {
        return params["orientation"].as<Annotation::Orientation>();
    }

    static LogLevel validateLogLevel(const po::variables_map& params) {
        const std::string logLevelStr = params["loglevel"].as<std::string>();

        if (logLevelStr == "DEBUG") {
            return LogLevel::DEBUG;
        } else if (logLevelStr == "INFO") {
            return LogLevel::INFO;
        } else if (logLevelStr == "WARNING") {
            return LogLevel::WARNING;
        } else if (logLevelStr == "ERROR") {
            return LogLevel::ERROR;
        } else {
            Logger::log(LogLevel::ERROR, "Invalid log level specified.");
            return LogLevel::INFO;
        }
    }
};
