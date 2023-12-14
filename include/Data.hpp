#ifndef DATA_HPP
#define DATA_HPP

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

#include "Align.hpp"
#include "Analysis.hpp"
#include "Clustering.hpp"
#include "Logger.hpp"
#include "SeqRickshaw.hpp"
#include "SplitReadCalling.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

using GroupsPath = std::map<std::string, fs::path>;  // path to ctrls, trtms

template <class DataType>
using GroupsMap = std::map<std::string, std::pair<std::string, DataType>>;

typedef std::vector<fs::path> PathVector;

class Data {
  /*
   * the type of the data:
   *
   */
 private:
  po::variables_map params;
  pt::ptree dataStructure;

  void createDirectories(fs::path &path);
  void createOutDir(fs::path &path, std::string &subcall);
  bool withControlData();

  GroupsPath retrieveGroupsPath(std::optional<fs::path> ctrls, fs::path trtms);
  void retrieveDataStructure(const GroupsPath &groupsPath);
  fs::path replaceParentDirPath(fs::path replacement, fs::path original);

 public:
  Data(po::variables_map _params);

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

  template <typename Callable>
  void callInAndOut(Callable f);

  // helper methods
  // determine content of directory and sort it (return as vector)
  PathVector getSortedDirContent(fs::path _path);
  // filter content of directory to only include files containing search string
  PathVector filterDirContent(PathVector vec, std::string sestr);
  std::string addSuffix(std::string _file, std::string _suffix, std::vector<std::string> _keys);

  // test stuff
  template <typename Callable>
  void bla(Callable f);

  void preproc();
  void align();
  void splitReadCalling();
  void clustering();
  void analysis();
};

#endif  // DATA_HPP
