#include "Data.hpp"

#include "Helper.hpp"

Data::Data(po::variables_map _params) : params(_params) {}

/**
 * Creates directories for the given path.
 *
 * @param path The path for which directories need to be created.
 */
void Data::createDirectories(fs::path &path) {
  if (fs::exists(path)) {
    Logger::log(LogLevel::INFO, "The direcotory", path.string(), " already exists.");
  } else {
    fs::create_directories(path);
    Logger::log(LogLevel::INFO, "Created directory ", path.string());
  }
}

/**
 * @brief Creates an output directory.
 *
 * This function creates an output directory at the specified path. The path is passed by reference
 * and is modified to contain the full path of the created directory. The subcall parameter is a
 * string that specifies additional information about the subcall being made.
 *
 * @param path The path where the output directory should be created.
 * @param subcall Additional information about the subcall being made.
 */
void Data::createOutDir(fs::path &path, std::string &subcall) {
  createDirectories(path);
  fs::path pathSub = path / subcall;
  createDirectories(pathSub);
}

/**
 * Replaces the parent directory path of a given file path with a new path.
 *
 * @param replacement The new parent directory path.
 * @param original The original file path.
 * @return The modified file path with the replaced parent directory path.
 */
fs::path Data::replaceParentDirPath(fs::path replacement, fs::path original) {
  return replacement / original.filename();
}

bool Data::withControlData() { return params["ctrls"].as<std::string>() != ""; }

pt::ptree Data::getDataStructure() { return this->dataStructure; }

void Data::setSubcall(std::string subcall) { params.at("subcall").value() = subcall; }

void Data::prepareSubcall(std::string subcall) {
  fs::path resultsDir = fs::path(params["outdir"].as<std::string>());
  createOutDir(resultsDir, subcall);

  if (subcall == "preproc") {
    preprocDataPrep();
  }

  if (subcall == "align") {
    alignDataPrep();
  }

  if (subcall == "detect") {
    detectDataPrep();
  }

  if (subcall == "clustering") {
    clusteringDataPrep();
  }

  if (subcall == "analysis") {
    analysisDataPrep();
  }
}

void Data::preprocDataPrep() {
  Logger::log(LogLevel::INFO, "Retrieving the data for preprocessing.");

  std::optional<fs::path> ctrlsPath = std::nullopt;

  if (withControlData()) {
    ctrlsPath = fs::path(this->params["ctrls"].as<std::string>());
  }
  fs::path trtmsPath = fs::path(this->params["trtms"].as<std::string>());

  GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);  //
  retrieveDataStructure(group);
}

void Data::alignDataPrep() {
  Logger::log(LogLevel::INFO, "Retrieving the data for alignment.");

  bool withPreprocessing = params["preproc"].as<std::bitset<1>>() == std::bitset<1>(1);

  std::optional<fs::path> ctrlsPath = std::nullopt;
  fs::path trtmsPath;

  if (withPreprocessing) {
    if (withControlData()) {
      ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "preproc/ctrls";
    }
    trtmsPath = fs::path(params["outdir"].as<std::string>()) / "preproc/trtms";
  } else {
    ctrlsPath = fs::path(params["ctrls"].as<std::string>());
    trtmsPath = fs::path(params["trtms"].as<std::string>());
  }

  GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
  retrieveDataStructure(group);
}

void Data::detectDataPrep() {
  Logger::log(LogLevel::INFO, "Retrieving the data for detecting split reads.");

  if (params["stats"].as<std::bitset<1>>() == 1) {
    fs::path statsfile = fs::path(params["outdir"].as<std::string>()) / "detect/detectStat.txt";
    std::ofstream ofs;
    ofs.open(statsfile.string());
    ofs << "sample\tmapped\tsplits\tmultisplits" << std::endl;
    ofs.close();
  }

  std::optional<fs::path> ctrlsPath = std::nullopt;
  if (withControlData()) {
    ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "align/ctrls";
  }
  fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "align/trtms";

  GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
  retrieveDataStructure(group);
}

void Data::clusteringDataPrep() {
  Logger::log(LogLevel::INFO, "Retrieving the data for clustering.");

  std::optional<fs::path> ctrlsPath = std::nullopt;

  if (withControlData()) {
    fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
  }
  fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

  GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
  retrieveDataStructure(group);
}

void Data::analysisDataPrep() {
  Logger::log(LogLevel::INFO, "Retrieving the data for analysis.");

  std::optional<fs::path> ctrlsPath = std::nullopt;

  if (withControlData()) {
    fs::path ctrlsPath = fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
  }

  fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

  GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
  retrieveDataStructure(group);
}

GroupsPath Data::retrieveGroupsPath(std::optional<fs::path> ctrls, fs::path trtms) {
  GroupsPath groups;

  if (!fs::is_directory(trtms)) {
    Logger::log(LogLevel::ERROR, "Not a directory: ", trtms.string());
    exit(EXIT_FAILURE);
  }

  groups.insert(std::make_pair("trtms", trtms));

  // Controls are optional
  if (!ctrls.has_value()) {
    Logger::log(LogLevel::WARNING, "This call runs RNAnue without control data");
    return groups;
  }

  if (!fs::is_directory(ctrls.value())) {
    Logger::log(LogLevel::ERROR, "Not a directory: ", ctrls.value().string());
    exit(EXIT_FAILURE);
  }

  groups.insert(std::make_pair("ctrls", ctrls.value()));

  return groups;
}

// creates property of the files in the system
void Data::retrieveDataStructure(const GroupsPath &groupsPath) {
  pt::ptree subcall;
  pt::ptree group;
  pt::ptree condition;

  for (const auto &groupPath : groupsPath) {
    Logger::log(LogLevel::INFO, "Retrieving data for group ", groupPath.first, " from ",
                groupPath.second.string());

    PathVector conditionsVec = getSortedDirContent(groupPath.second);

    for (const auto &conditionPath : conditionsVec) {
      if (!fs::is_directory(conditionPath)) {
        Logger::log(LogLevel::WARNING,
                    "Not a directory (will be ignored): ", conditionPath.string());
        continue;
      }
      condition = retrieveConditionTree(groupPath.first, conditionPath);

      group.push_back(std::make_pair("", condition));
    }
    subcall.add_child(groupPath.first, group);
    group.clear();
  }

  dataStructure.add_child(params["subcall"].as<std::string>(), subcall);

  fs::path dataStructurePath = params["outdir"].as<std::string>() / fs::path("data.json");

  Logger::log(LogLevel::INFO, "Saving parsed data structure to file ", dataStructurePath.string(),
              ".");
  boost::property_tree::json_parser::write_json(dataStructurePath.string(), dataStructure);
}

pt::ptree Data::retrieveConditionTree(std::string group, fs::path conditionPath) {
  std::string subcall = params["subcall"].as<std::string>();

  int expectedElementCount = 0;
  int elementCounter = 0;

  pt::ptree condition, samples;

  std::vector<std::string> sampleKeys;
  PathVector dataFiles = getSortedDirContent(conditionPath);

  std::string conditionName = conditionPath.stem().string();
  Logger::log(LogLevel::INFO, "Retrieving data for condition ", conditionName, " from ",
              conditionPath.string());

  if (subcall == "preproc") {
    expectedElementCount = (params["readtype"].as<std::string>() == "PE") ? 2 : 1;
    sampleKeys = {"forward", "reverse"};
    dataFiles = filterDirContent(dataFiles, ".fastq");
  } else if (subcall == "align") {
    expectedElementCount = 1;
    sampleKeys = {"forward"};
    if (params["preproc"].as<std::bitset<1>>() == std::bitset<1>(1)) {
      dataFiles = filterDirContent(dataFiles, "_preproc.fastq");
    } else {
      dataFiles = filterDirContent(dataFiles, ".fastq");
    }
  } else if (subcall == "detect") {
    expectedElementCount = 1;
    sampleKeys = {"matched"};
    dataFiles = filterDirContent(dataFiles, "matched.sam");
  } else if (subcall == "clustering" || subcall == "analysis") {
    expectedElementCount = 1;
    sampleKeys = {"splits"};
    dataFiles = filterDirContent(dataFiles, "_splits.sam");
  }

  fs::path outConditionDir = fs::path(params["outdir"].as<std::string>()) /
                             params["subcall"].as<std::string>() / group / conditionPath.filename();

  pt::ptree sample, files;
  for (const auto &dataFile : dataFiles) {
    Logger::log(LogLevel::INFO, "Retrieving data for sample ", dataFile.string());

    files.put(sampleKeys[elementCounter % expectedElementCount], dataFile.string());

    if (elementCounter == expectedElementCount - 1) {
      sample.add_child("input", files);

      pt::ptree output = retrieveSampleOutputTree(outConditionDir, sample);
      sample.add_child("output", output);

      samples.push_back(std::make_pair("", sample));
      elementCounter = 0;
      files.clear();
      sample.clear();
    } else {
      ++elementCounter;
    }
  }

  condition.put("condition", conditionName);
  condition.add_child("samples", samples);

  return condition;
}

void printTree(const boost::property_tree::ptree &pt, int level = 0) {
  for (const auto &node : pt) {
    std::cout << std::string(level * 2, ' ') << node.first << ": "
              << node.second.get_value<std::string>() << "\n";
    printTree(node.second, level + 1);
  }
}

// create ptree of the output files
pt::ptree Data::retrieveSampleOutputTree(fs::path outConditionDir, pt::ptree inputTree) {
  pt::ptree output;

  if (params["subcall"].as<std::string>() == "preproc") {
    // replace input path with output path (results/...)
    printTree(inputTree);

    fs::path fwd = fs::path(inputTree.get<std::string>("input.forward"));
    std::string forward = replaceParentDirPath(outConditionDir, fwd).string();

    // create output files using the input files with an additonal suffix
    std::string forwardOut = addSuffix(forward, "_preproc", {"_R1", "fwd"});
    output.put("forward", forwardOut);  // write output back to ptree

    // using paired-end reads results in additional files for unassembled
    // reads
    if (params["readtype"].as<std::string>() == "PE") {
      // replace
      fs::path rev = fs::path(inputTree.get<std::string>("input.reverse"));
      std::string reverse = replaceParentDirPath(outConditionDir, rev).string();

      // create outfiles using the input files with an additional suffix
      std::string r1only = addSuffix(reverse, "_R1only", {"_R2", "rev"});
      std::string r2only = addSuffix(reverse, "_R2only", {"_R2", "rev"});

      std::string r1unmerged = addSuffix(reverse, "_R1unmerged", {"_R2", "rev"});
      std::string r2unmerged = addSuffix(reverse, "_R2unmerged", {"_R2", "rev"});

      output.put("R1only", r1only);  // push both back to output ptree
      output.put("R2only", r2only);

      output.put("R1unmerged", r1unmerged);
      output.put("R2unmerged", r2unmerged);
    }
  } else if (params["subcall"].as<std::string>() == "align") {
    fs::path fwdInPath =
        fs::path(inputTree.get<std::string>("input.forward")).replace_extension(".sam");
    std::string fwdOutPath = replaceParentDirPath(outConditionDir, fwdInPath).string();
    std::string matched = addSuffix(fwdOutPath, "_matched", {});

    output.put("matched", matched);
  } else if (params["subcall"].as<std::string>() == "detect") {
    fs::path matchedInPath = fs::path(inputTree.get<std::string>("input.matched"));
    std::string splitsGeneralPath = replaceParentDirPath(outConditionDir, matchedInPath).string();
    std::string splitsOutPath = addSuffix(splitsGeneralPath, "splits", {"matched"});
    std::string multsplitsOutPath = addSuffix(splitsGeneralPath, "multsplits", {"matched"});

    output.put("splits", splitsOutPath);
    output.put("multsplits", multsplitsOutPath);
  } else if (params["subcall"].as<std::string>() == "clustering") {
    fs::path splits = fs::path(inputTree.get<std::string>("input.splits"));
    splits.replace_extension(".txt");
    std::string clu = replaceParentDirPath(outConditionDir, splits).string();
    std::string clusters = addSuffix(clu, "_clusters", {"_splits"});
    // output.put("clusters", clusters); # no individual outputs for each -
    // instead one file
  } else if (params["subcall"].as<std::string>() == "analysis") {
    fs::path splits = fs::path(inputTree.get<std::string>("input.splits"));
    splits.replace_extension(".txt");
    std::string its = replaceParentDirPath(outConditionDir, splits).string();
    std::string ints = addSuffix(its, "_interactions", {"_splits"});
    output.put("interactions", ints);
  }
  return output;
}

//
template <typename Callable>
void Data::callInAndOut(Callable f) {
  // retrieve paths for parameters
  std::string subcallStr = params["subcall"].as<std::string>();
  prepareSubcall(subcallStr);

  fs::path outDir = fs::path(params["outdir"].as<std::string>());
  fs::path outSubcallDir = outDir / fs::path(subcallStr);

  std::cout << "Calling " << subcallStr << " start creating output dir" << std::endl;

  // create output directory (and subdirectory for subcall)
  createDirectories(outDir);         // create output
  createDirectories(outSubcallDir);  //  create subcall

  pt::ptree subcall;

  try {
    subcall = dataStructure.get_child(subcallStr);
  } catch (pt::ptree_error &e) {
    std::cout << "### ERROR - " << subcallStr << " has not been found in the data structure"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::deque<std::string> groups = {"trtms"};
  if (subcall.size() > 1) {
    groups.push_front("ctrls");
  }

  pt::ptree conditions, samples, path;
  fs::path outGroupDir, outConditionDir;

  std::cout << "*** create directories to store the results (specified via --outdir)" << std::endl;

  for (unsigned i = 0; i < groups.size(); ++i) {
    // create directory for groups (e.g., ctrls, trtms)
    outGroupDir = outSubcallDir / fs::path(groups[i]);

    // don't create subfolders for clustering && analysis
    // if(params["subcall"].as<std::string>() != "clustering") {
    createDirectories(outGroupDir);
    //}

    conditions = subcall.get_child(groups[i]);
    BOOST_FOREACH (pt::ptree::value_type const &v, conditions.get_child("")) {
      pt::ptree condition = v.second;

      // create directory for condition (e.g., rpl_exp)
      outConditionDir = outGroupDir / fs::path(condition.get<std::string>("condition"));
      createDirectories(outConditionDir);

      samples = condition.get_child("samples");
      // iterate over samples
      BOOST_FOREACH (pt::ptree::value_type const &w, samples.get_child("")) {
        pt::ptree sample = w.second;

        // call start function
        f(sample);
        //
        // f(path, group  sample)
      }
    }
  }
}

template <typename Callable>
void Data::bla(Callable f) {
  std::cout << "bla test" << std::endl;
  f();
}

void Data::preproc() {
  SeqRickshaw srs(params);
  callInAndOut(std::bind(&SeqRickshaw::start, srs, std::placeholders::_1));
}

void Data::align() {
  Align aln(params);
  callInAndOut(std::bind(&Align::start, aln, std::placeholders::_1));
}

void Data::splitReadCalling() {
  SplitReadCalling src(params);
  callInAndOut(std::bind(&SplitReadCalling::start, src, std::placeholders::_1));
}

void Data::clustering() {
  // create Object of CLustering
  Clustering clu(params);
  callInAndOut(std::bind(&Clustering::start, &clu, std::placeholders::_1));

  clu.sumup();
}

void Data::analysis() {
  Analysis anl(params);
  callInAndOut(std::bind(&Analysis::start, &anl, std::placeholders::_1));

  if (params["outcnt"].as<std::bitset<1>>() == std::bitset<1>(1)) {
    anl.createCountTable();
  }
}

/*
 * Helper Function
 */

// lists and sorts the content of a directory
PathVector Data::getSortedDirContent(fs::path _path) {
  PathVector content;  //
  copy(fs::directory_iterator(_path), fs::directory_iterator(), back_inserter(content));
  sort(content.begin(), content.end());  // sort the content
  return content;
}

// filters content of directory
PathVector Data::filterDirContent(PathVector vec, std::string sestr) {
  PathVector content;  // create new vector
  std::copy_if(vec.begin(), vec.end(), std::back_inserter(content), [&sestr](fs::path x) {
    // only return if fs::path contains sestr
    return x.string().find(sestr) != std::string::npos;
  });
  return content;
}

// adds suffix to filename (optional:
std::string Data::addSuffix(std::string _file, std::string _suffix,
                            std::vector<std::string> _keywords) {
  int keyPos, tmpKeyPos = -1;    // buffer the positions of the keywords
  int dotPos = _file.find(".");  // determine position of the dot
  keyPos = dotPos;
  if (!_keywords.empty()) {
    for (unsigned i = 0; i < _keywords.size(); ++i) {
      tmpKeyPos = _file.find(_keywords[i]);
      if (tmpKeyPos != -1) {  // key could be found
        keyPos = tmpKeyPos;
      }
    }
  }
  std::string newFile = _file.substr(0, keyPos);
  newFile += _suffix;
  newFile += _file.substr(dotPos, _file.size());

  return newFile;
}