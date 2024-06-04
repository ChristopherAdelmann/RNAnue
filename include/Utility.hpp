#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan3
#include <seqan3/io/sam_file/all.hpp>

// Standard
#include <chrono>
#include <iomanip>
#include <iostream>
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
