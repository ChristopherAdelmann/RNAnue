//
// Created by Richard Albin Schaefer on 1/23/24.
//

#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <seqan3/io/sequence_file/all.hpp>
#include <unordered_set>

#include "DataTypes.hpp"

namespace fs = boost::filesystem;

// filesystem manipulation
namespace helper {
void createDir(fs::path path, std::ostream& out);
void createOutDir(fs::path path, std::ostream& out);
void createTmpDir(fs::path subpath);

std::size_t countUniqueSamEntries(fs::path path);
std::size_t countSamEntries(fs::path path);

dtp::PathVector listDirFiles(fs::path& path);
dtp::PathVector filterDirFiles(dtp::PathVector& pathvec, std::string subcall);
bool caseInsensitivePathCompare(const fs::path& a, const fs::path& b);
void renameFiles(fs::path dir, std::string extension);
dtp::PathVector genOutPath(dtp::PathVector pathvec, std::string dir);
void printTree(const boost::property_tree::ptree& pt, int level);

fs::path replacePath(fs::path _replacement, fs::path _original);
void deleteDir(fs::path path);
std::string addSuffix(std::string filename, std::string suffix, std::vector<std::string> keywords);

void splitFileN(fs::path source, fs::path targetPath, int sep);
int lineCount(fs::path file, std::string command);  // count lines of file
std::vector<fs::path> listFiles(fs::path path);
void concatFiles(fs::path target, std::vector<fs::path> files);

std::string getTime();  // reports the current time
void simulateProcessing(std::chrono::microseconds desiredDuration);
class Timer {
   public:
    Timer();
    ~Timer();
    void stop();

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};
}  // namespace helper

namespace seqIO {
bool filterReads(auto& qual, int quality, int minlen);

dtp::DNAVector spanToVector(
    dtp::DNASpan span);  // convert a dna5 span to a vector (for easier manipulation)
// determine the alphabet of a sequence (required to preprocess the search pattern)
dtp::DNAVector determineAlphabet(dtp::DNAVector seq);

// functions for debugging
void printStates(dtp::DNAVector seq, std::string state, dtp::Left left, dtp::Right right,
                 int readPos);
void printDNAVector(dtp::DNAVector seq,
                    std::ostream& ofs);  // prints a DNAVector to a stream (e.g., file)p
void printDNASpan(dtp::DNASpan span, std::ostream& ofs);
void printQualSpan(dtp::QualSpan span, std::ostream& ofs);

}  // namespace seqIO

#endif  // UTILITY_HPP
