//
// Created by Richard Albin Schaefer on 1/24/24.
//

#ifndef RNANUE_DATA_HPP
#define RNANUE_DATA_HPP

#include <math.h>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <typeinfo>

#include "Align.hpp"
#include "Analysis.hpp"
#include "Clustering.hpp"
#include "Constansts.hpp"
#include "Logger.hpp"
#include "Preprocess.hpp"
#include "SplitReadCalling.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;
namespace jp = boost::property_tree::json_parser;
namespace pi = constants::pipelines;

// map that contains the paths to the groups (e.g., ctrls, trtms)
using GroupsPath = std::map<std::string, fs::path>;
using PathVector = std::vector<fs::path>;

class Data {
   public:
    Data(po::variables_map params);
    ~Data() = default;
    //
    template <typename Callable>
    void callInAndOut(Callable f);

   private:
    po::variables_map params;
    pt::ptree dataStructure;

    void createDirectories(fs::path& path);
    void createOutDir(fs::path& path, std::string& subcall);
    bool withControlData();

    GroupsPath retrieveGroupsPath(std::optional<fs::path> ctrls, fs::path trtms);
    void retrieveDataStructure(const GroupsPath& groupsPath);
    fs::path replaceParentDirPath(fs::path replacement, fs::path original);

   public:
    // getter & setter
    //
    pt::ptree getDataStructure();
    void setSubcall(std::string subcall);
    void prepareSubcall(std::string subcall);

    //
    void preprocDataPrep();
    void alignDataPrep();
    void detectDataPrep();
    void clusteringDataPrep();
    void analysisDataPrep();

    //
    pt::ptree retrieveConditionTree(std::string group, fs::path _condition);
    pt::ptree retrieveSampleOutputTree(fs::path outConditionDir, pt::ptree inputTree);

    // iterates through the data structure and excutes the subcall

    // helper methods
    // determine content of directory and sort it (return as vector)
    PathVector getSortedDirContent(fs::path _path);
    // filter content of directory to only include files containing search string
    PathVector filterDirContent(PathVector vec, std::string sestr);
    std::string addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keys);

    void preproc();
    void align();
    void splitReadCalling();
    void clustering();
    void analysis();
};

#endif  // RNANUE_DATA_HPP
