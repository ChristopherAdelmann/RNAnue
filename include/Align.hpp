#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// Standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <optional>
#include <regex>
#include <sstream>
#include <utility>

// Class
#include "Logger.hpp"
#include "Utility.hpp"

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

class Align {
   public:
    explicit Align(po::variables_map params);
    ~Align() = default;

    // alignment
    void start(pt::ptree sample);

   private:
    po::variables_map params;
    std::string index;
    std::string segemehlSysCall;

    void sortAlignments(const std::string &alignmentsPath);
    void buildIndex();
    void alignReads(const std::string &query, const std::string &mate, const std::string &matched);
};
