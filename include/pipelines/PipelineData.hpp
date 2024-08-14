#pragma once

// Standard
#include <array>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

// Classes
#include "Logger.hpp"

namespace pipelines {

namespace fs = std::filesystem;

template <typename Container>
concept StringContainer = std::ranges::range<Container> &&
                          std::same_as<std::ranges::range_value_t<Container>, std::string>;

struct PipelineData {
   protected:
    static bool isHidden(const std::filesystem::path &path) {
        std::string filename = path.filename().string();
        return !filename.empty() && filename[0] == '.';
    }

    static std::string getSampleName(const fs::path &filePath) {
        std::string fileName = filePath.stem().string();
        size_t underscorePos = fileName.find("_");
        if (underscorePos != std::string::npos) {
            return fileName.substr(0, underscorePos);
        } else {
            return fileName;
        }
    }

    static bool hasSuffix(const std::string &fullString, const std::string &ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(),
                                            ending));
        } else {
            return false;
        }
    };

    template <StringContainer Container>
    static bool hasAnySuffix(const std::string &fullString, const Container &endings) {
        for (const auto &ending : endings) {
            if (hasSuffix(fullString, ending)) {
                return true;
            }
        }
        return false;
    };

    static bool hasPrefix(const std::string &fullString, const std::string &prefix) {
        if (fullString.length() >= prefix.length()) {
            return (0 == fullString.compare(0, prefix.length(), prefix));
        } else {
            return false;
        }
    };

    template <StringContainer Container>
    static bool hasAnyPrefix(const std::string &fullString, const Container &prefixes) {
        for (const auto &prefix : prefixes) {
            if (hasPrefix(fullString, prefix)) {
                return true;
            }
        }
        return false;
    };

    static bool contains(const std::string &fullString, const std::string &substring) {
        return fullString.find(substring) != std::string::npos;
    };

    template <StringContainer Container>
    static bool containsAny(const std::string &fullString, const Container &substrings) {
        for (const auto &substring : substrings) {
            if (contains(fullString, substring)) {
                return true;
            }
        }
        return false;
    };

    static std::vector<fs::path> getDirectories(const fs::path &parentDir) {
        std::vector<fs::path> directories;
        for (const auto &entry : fs::directory_iterator(parentDir)) {
            if (entry.is_directory()) {
                directories.push_back(entry.path());
            } else if (entry.is_regular_file() && !isHidden(entry)) {
                Logger::log(LogLevel::WARNING,
                            "Found file in directory, which is not allowed: ", entry.path());
            }
        }
        return directories;
    };

    template <StringContainer Container = std::vector<std::string>>
    static std::vector<fs::path> getValidFilePaths(const fs::path &directory,
                                                   const Container &allowedSuffixes = Container{},
                                                   const Container &allowedPrefixes = Container{}) {
        std::vector<fs::path> filePaths;

        for (const auto &entry : fs::directory_iterator(directory)) {
            if (!entry.is_regular_file()) {
                Logger::log(LogLevel::WARNING,
                            "Found not supported type in directory: ", entry.path());
                continue;
            }

            const auto &filePathStr = entry.path().string();

            bool hasValidSuffix = hasAnySuffix(filePathStr, allowedSuffixes);
            bool hasValidPrefix = hasAnyPrefix(filePathStr, allowedPrefixes);

            if (hasValidSuffix && hasValidPrefix) {
                filePaths.push_back(entry.path());
            }
        }

        return filePaths;
    };
};
}  // namespace pipelines
