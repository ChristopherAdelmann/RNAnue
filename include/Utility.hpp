#pragma once

// Standard
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

// Boost
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

// Internal
#include "Logger.hpp"

namespace helper {

constexpr int DOUBLE_COMPARISON_GRACE_FACTOR = 10;

inline auto isEqual(double lhs, double rhs,
                    double epsilon = std::numeric_limits<double>::epsilon() *
                                     DOUBLE_COMPARISON_GRACE_FACTOR) -> bool {
    return fabs(lhs - rhs) < epsilon;
}

template <typename Container>
concept StringContainer = std::ranges::range<Container> &&
                          std::same_as<std::ranges::range_value_t<Container>, std::string>;

namespace fs = std::filesystem;

void createTmpDir(const fs::path &subpath);
void deleteDir(const fs::path &path);

inline auto hasSuffix(const std::string &fullString, const std::string &ending) -> bool {
    if (fullString.length() >= ending.length()) {
        return (0 ==
                fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    return false;
};
template <StringContainer Container>
auto hasAnySuffix(const std::string &fullString, const Container &endings) -> bool {
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

inline auto hasPrefix(const std::string &fullString, const std::string &prefix) -> bool {
    if (fullString.length() >= prefix.length()) {
        return (0 == fullString.compare(0, prefix.length(), prefix));
    }
    return false;
};

template <StringContainer Container>
auto hasAnyPrefix(const std::string &fullString, const Container &prefixes) -> bool {
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
auto getValidFilePaths(const fs::path &directory,
                       const Container &allowedSuffixes = Container{},  // NOLINT
                       const Container &allowedPrefixes = Container{}) -> std::vector<fs::path> {
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

auto getDirIfExists(const fs::path &path) -> std::optional<fs::path>;

auto getUUID() -> std::string;

void mergeSamFiles(const std::vector<fs::path> &inputPaths, const fs::path &outputPath);
void mergeFastqFiles(const std::vector<fs::path> &inputPaths, const fs::path &outputPath);

auto getFilePathsInDir(const fs::path &dir) -> std::vector<fs::path>;

void concatAndDeleteFilesInTmpDir(const fs::path &tmpDir, const fs::path &outPath);

auto countUniqueSamEntries(fs::path &path) -> std::size_t;
auto countSamEntries(fs::path &path) -> std::size_t;
auto countSamEntriesSeqAn(fs::path &path) -> std::size_t;

/** Checks if a value is contained within a specified range, with a given tolerance.
 *
 * @param value The value to check.
 * @param comparisonValue The value to compare against.
 * @param tolerance The tolerance within which the values are considered equal.
 * @return true if the value is contained within the range, false otherwise.
 **/
auto isContained(int32_t value, int32_t comparisonValue, int32_t tolerance) -> bool;

template <typename T>
auto calculateMedian(std::vector<T> values) -> T {
    std::sort(values.begin(), values.end());
    const auto size = values.size();
    if (size % 2 == 0) {
        return (values[size / 2 - 1] + values[size / 2]) / 2;
    } else {
        return values[size / 2];
    }
}

auto generateRandomHexColor() -> std::string;

auto getTime() -> std::string;  // reports the current time

class Timer {
   public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    Timer(const Timer &) = default;
    Timer(Timer &&) = delete;
    auto operator=(const Timer &) -> Timer & = default;
    auto operator=(Timer &&) -> Timer & = delete;
    ~Timer() { stop(); }

    void stop();

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};
}  // namespace helper
