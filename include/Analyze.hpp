#pragma once

// Standard
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <vector>

// Boost
// #include <boost/accumulators/accumulators.hpp>
// #include <boost/accumulators/statistics.hpp>
#include <boost/filesystem.hpp>
// #include <boost/math/distributions/binomial.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// SeqAn3
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>

// Class
// #include "IBPTree.hpp"

// define tags
using seqan3::operator""_tag;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;
// namespace ma = boost::math;

class Analyze {
   public:
    explicit Analyze(po::variables_map params);
    ~Analyze() = default;

    void start(pt::ptree sample);
    // void processSingle(std::string& single);
    // void processSplits(std::string& splits, std::string& output);
    // void normalize();  // normalize the frequencies to 1

    // // write output files (of the analysis)
    // void writeInteractionsHeader(std::ofstream& fout);
    // void writeAllIntsHeader(std::ofstream& fout);
    // void addToAllIntsHeader(std::ofstream& fout, std::string key);
    // void writeAllInts();

    // // other operations
    // void addToFreqMap(std::pair<std::string, std::string> key, double value);
    // void addToSuppReads(dtp::IntPartner p1, dtp::IntPartner p2);
    // void addToCompl(dtp::IntPartner p1, dtp::IntPartner p2, double complementarity);
    // void addToHybEnergies(dtp::IntPartner p1, dtp::IntPartner p2, double hybenergy);

    // // calculate filters
    // double calcGCS(std::vector<double>& complementarities);
    // double calcGHS(std::vector<double>& hybenergies);
    // double calcStat(dtp::IntKey key, int x);

   private:
    po::variables_map params;
    // IBPTree features;
    // std::map<std::pair<std::string, std::string>, double> freq;  // strand, name
    // std::string condition;                                       // buffers the current condition
    // // maps for storing filters and suppreads
    // std::map<dtp::IntKey, std::vector<double>> suppreads;
    // std::map<dtp::IntKey, std::vector<std::vector<double>>> complementarities;
    // std::map<dtp::IntKey, std::vector<std::vector<double>>> hybenergies;

    int repcount;   // counter for current replicate
    int readcount;  // total number of reads
};