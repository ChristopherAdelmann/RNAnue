#pragma once

// Standard
#include <concepts>
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <filesystem>
#include <string>

// Classes
#include "Logger.hpp"

namespace po = boost::program_options;

template <typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

struct ParameterValidator {
    template <typename T>
        requires arithmetic<T>
    static T validateArithmetic(const po::variables_map& params, const std::string& paramName,
                                const T lowerLimit, const T upperLimit) {
        const int value = params[paramName].as<T>();

        if (value <= lowerLimit && value >= upperLimit) {
            Logger::log(LogLevel::ERROR, paramName + " must be an integer between " +
                                             std::to_string(lowerLimit) + " and " +
                                             std::to_string(upperLimit));
        }

        return value;
    }

    static std::filesystem::path validateFilePath(const po::variables_map& params,
                                                  const std::string& paramName) {
        const std::string filePathStr = params[paramName].as<std::string>();
        std::filesystem::path filePath = std::filesystem::path(filePathStr);

        if (!std::filesystem::exists(filePath) || std::filesystem::is_directory(filePath)) {
            Logger::log(LogLevel::ERROR, paramName + " is not a valid file path.");
        }

        return filePath;
    }

    static std::filesystem::path validateDirectory(const po::variables_map& params,
                                                   const std::string& paramName) {
        const std::string dirPathStr = params[paramName].as<std::string>();
        std::filesystem::path dirPath = std::filesystem::path(dirPathStr);

        if (!std::filesystem::is_directory(dirPath)) {
            Logger::log(LogLevel::ERROR, paramName + " is not a valid directory.");
        }

        return dirPath;
    }
};
