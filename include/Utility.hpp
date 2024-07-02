#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>

// Standard
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <unordered_set>

// Classes
#include "Constants.hpp"
#include "DataTypes.hpp"
#include "Logger.hpp"

namespace pi = constants::pipelines;
namespace fs = boost::filesystem;

// filesystem manipulation
namespace helper {
void createTmpDir(fs::path subpath);
void deleteDir(fs::path path);

void mergeSamFiles(std::vector<fs::path> inputPaths, fs::path outputPath);

std::size_t countUniqueSamEntries(fs::path path);
std::size_t countSamEntries(fs::path path);
std::size_t countSamEntriesSeqAn(fs::path path);

void printTree(const boost::property_tree::ptree& pt, int level);

template <typename T>
T calculateMedian(std::vector<T> values) {
    std::sort(values.begin(), values.end());
    const auto size = values.size();
    if (size % 2 == 0) {
        return (values[size / 2 - 1] + values[size / 2]) / 2;
    } else {
        return values[size / 2];
    }
};

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
