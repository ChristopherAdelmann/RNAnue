#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

// Standard
#include <deque>
#include <iostream>
#include <optional>

// Classes
#include "Align.hpp"
#include "Analyze.hpp"
#include "Cluster.hpp"
#include "Constants.hpp"
#include "Detect.hpp"
#include "Logger.hpp"
#include "Preprocess.hpp"

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
    explicit Data(po::variables_map params);
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
    template <typename... Seeds>
    PathVector filterDirContent(const PathVector& paths, Seeds... seeds);
    std::string addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keys);

    void preproc();
    void align();
    void splitReadCalling();
    void clustering();
    void analysis();
};