#pragma once

// Standard
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// Boost
#include <boost/property_tree/ptree.hpp>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>

// Classes
#include "Constants.hpp"
#include "DataTypes.hpp"
#include "Logger.hpp"

namespace helper {

template <typename Container>
concept StringContainer = std::ranges::range<Container> &&
                          std::same_as<std::ranges::range_value_t<Container>, std::string>;

namespace fs = std::filesystem;

void createTmpDir(fs::path subpath);
void deleteDir(fs::path path);

inline bool hasSuffix(const std::string &fullString, const std::string &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 ==
                fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
};
template <StringContainer Container>
bool hasAnySuffix(const std::string &fullString, const Container &endings) {
    if (endings.empty()) {
        return true;
    }

    for (const auto &ending : endings) {
        if (hasSuffix(fullString, ending)) {
            return true;
        }
    }
    return false;
};

inline bool hasPrefix(const std::string &fullString, const std::string &prefix) {
    if (fullString.length() >= prefix.length()) {
        return (0 == fullString.compare(0, prefix.length(), prefix));
    } else {
        return false;
    }
};

template <StringContainer Container>
bool hasAnyPrefix(const std::string &fullString, const Container &prefixes) {
    if (prefixes.empty()) {
        return true;
    }

    for (const auto &prefix : prefixes) {
        if (hasPrefix(fullString, prefix)) {
            return true;
        }
    }
    return false;
};

template <StringContainer Container = std::vector<std::string>>
std::vector<fs::path> getValidFilePaths(const fs::path &directory,
                                        const Container &allowedSuffixes = Container{},
                                        const Container &allowedPrefixes = Container{}) {
    std::vector<fs::path> filePaths;

    for (const auto &entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file()) {
            Logger::log(LogLevel::WARNING, "Found not supported type in directory: ", entry.path());
            continue;
        }

        Logger::log(LogLevel::DEBUG, "Found file: ", entry);

        const auto &filePathStr = entry.path().string();

        bool hasValidSuffix = hasAnySuffix(filePathStr, allowedSuffixes);
        bool hasValidPrefix = hasAnyPrefix(filePathStr, allowedPrefixes);

        if (hasValidSuffix && hasValidPrefix) {
            filePaths.push_back(entry.path());
        }
    }

    return filePaths;
};

std::optional<fs::path> getDirIfExists(const fs::path &path);

void mergeSamFiles(std::vector<fs::path> inputPaths, fs::path outputPath);

std::size_t countUniqueSamEntries(fs::path path);
std::size_t countSamEntries(fs::path path);
std::size_t countSamEntriesSeqAn(fs::path path);

/** Checks if a value is contained within a specified range, with a given tolerance.
 *
 * @param value The value to check.
 * @param comparisonValue The value to compare against.
 * @param tolerance The tolerance within which the values are considered equal.
 * @return true if the value is contained within the range, false otherwise.
 **/
bool isContained(const int32_t value, const int32_t comparisonValue, const int32_t tolerance);

void printTree(const boost::property_tree::ptree &pt, int level);

template <typename T>
T calculateMedian(std::vector<T> values) {
    std::sort(values.begin(), values.end());
    const auto size = values.size();
    if (size % 2 == 0) {
        return (values[size / 2 - 1] + values[size / 2]) / 2;
    } else {
        return values[size / 2];
    }
}

std::string generateRandomHexColor();

std::string getTime();  // reports the current time
class Timer {
   public:
    Timer();
    ~Timer();
    void stop();

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};
}  // namespace helper
