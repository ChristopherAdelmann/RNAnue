#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <math.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <fstream>
#include <iostream>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <string>
#include <tuple>
#include <vector>

#include "CustomSamTags.hpp"
#include "Logger.hpp"

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef std::pair<std::pair<int, int>, std::string> Feature;

class Analyze {
   private:
    po::variables_map params;
    std::map<std::string, std::vector<std::pair<std::pair<int, int>, std::string>>> features;

    std::vector<std::string> interPaths;  // paths of the interactions file

   public:
    Analyze(po::variables_map params);

    std::string retrieveTagValue(std::string tags, std::string tagName, std::string oldValue);

    void createCountTable();

    void parseAnnotations();
    void start(pt::ptree sample);
};

#endif  // ANALYSIS_HPP
#ifndef RNANUE_ANALYSIS_HPP
#define RNANUE_ANALYSIS_HPP

#include <iostream>

#endif  // RNANUE_ANALYSIS_HPP
