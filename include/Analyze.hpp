#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// Standard
#include <math.h>

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/utility/views/chunk.hpp>

// Classes
#include "CustomSamTags.hpp"
#include "FeatureAnnotator.hpp"
#include "Logger.hpp"

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef std::pair<std::pair<int, int>, std::string> Feature;

class Analyze {
   public:
    explicit Analyze(po::variables_map params);
    ~Analyze() = default;

    void start(pt::ptree sample);

    void createCountTable();

   private:
    po::variables_map params;
    std::map<std::string, std::vector<std::pair<std::pair<int, int>, std::string>>> features;
    Annotation::FeatureAnnotator annotator;

    std::vector<std::string> interPaths;  // paths of the interactions file

    std::string retrieveTagValue(std::string tags, std::string tagName, std::string oldValue);
};
