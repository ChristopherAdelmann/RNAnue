#include "Data.hpp"

#include "Utility.hpp"

Data::Data(po::variables_map _params) : params(_params) {}

/**
 * Creates directories for the given path.
 *
 * @param path The path for which directories need to be created.
 */
void Data::createDirectories(const fs::path &path) {
    if (fs::exists(path)) {
        Logger::log(LogLevel::INFO, "The directory ", path.string(), " already exists");
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
void Data::createOutDir(const fs::path &path, const std::string &subcall) {
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

bool Data::withControlData() { return !params["ctrls"].as<std::string>().empty(); }

pt::ptree Data::getDataStructure() { return this->dataStructure; }

void Data::setSubcall(std::string subcall) { params.at("subcall").value() = subcall; }

void Data::prepareSubcall(std::string subcall) {
    fs::path resultsDir = fs::path(params["outdir"].as<std::string>());
    createOutDir(resultsDir, subcall);

    if (subcall == pi::PREPROCESS) {
        preprocDataPrep();
    }

    if (subcall == pi::ALIGN) {
        alignDataPrep();
    }

    if (subcall == pi::DETECT) {
        detectDataPrep();
    }

    if (subcall == pi::ANALYZE) {
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

    bool withPreprocessing = params[pi::PREPROCESS].as<bool>();

    std::optional<fs::path> ctrlsPath = std::nullopt;
    fs::path trtmsPath;

    if (withPreprocessing) {
        if (withControlData()) {
            ctrlsPath = fs::path(params["outdir"].as<std::string>()) / pi::PREPROCESS / "ctrls";
        }
        trtmsPath = fs::path(params["outdir"].as<std::string>()) / pi::PREPROCESS / "trtms";
    } else {
        ctrlsPath = fs::path(params["ctrls"].as<std::string>());
        trtmsPath = fs::path(params["trtms"].as<std::string>());
    }

    GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
    retrieveDataStructure(group);
}

void Data::detectDataPrep() {
    Logger::log(LogLevel::INFO, "Retrieving the data for detecting split reads.");

    std::optional<fs::path> ctrlsPath = std::nullopt;
    if (withControlData()) {
        ctrlsPath = fs::path(params["outdir"].as<std::string>()) / pi::ALIGN / "ctrls";
    }
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / pi::ALIGN / "trtms";

    GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
    retrieveDataStructure(group);
}

void Data::analysisDataPrep() {
    Logger::log(LogLevel::INFO, "Retrieving the data for analysis.");

    std::optional<fs::path> ctrlsPath = std::nullopt;

    if (withControlData()) {
        fs::path(params["outdir"].as<std::string>()) / "detect/ctrls";
    }
    fs::path trtmsPath = fs::path(params["outdir"].as<std::string>()) / "detect/trtms";

    GroupsPath group = retrieveGroupsPath(ctrlsPath, trtmsPath);
    retrieveDataStructure(group);
}

GroupsPath Data::retrieveGroupsPath(std::optional<fs::path> ctrls, fs::path trtms) {
    GroupsPath groups;

    if (!fs::is_directory(trtms)) {
        Logger::log(LogLevel::ERROR, "Not a directory: ", trtms.string());
    }

    groups.insert(std::make_pair("trtms", trtms));

    // Controls are optional
    if (!ctrls.has_value()) {
        Logger::log(LogLevel::WARNING, "This call runs RNAnue without control data");
        return groups;
    }

    if (!fs::is_directory(ctrls.value())) {
        Logger::log(LogLevel::ERROR, "Not a directory: ", ctrls.value().string());
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

    size_t expectedElementCount = 0;
    size_t elementCounter = 0;

    pt::ptree condition, samples;

    std::vector<std::string> sampleKeys;
    PathVector dataFiles = getSortedDirContent(conditionPath);

    std::string conditionName = conditionPath.stem().string();
    Logger::log(LogLevel::INFO, "Retrieving data for condition ", conditionName, " from ",
                conditionPath.string());

    if (subcall == pi::PREPROCESS) {
        expectedElementCount = (params["readtype"].as<std::string>() == "PE") ? 2 : 1;
        sampleKeys = {"forward", "reverse"};
        dataFiles = filterDirContent(dataFiles, ".fastq");
    } else if (subcall == pi::ALIGN) {
        expectedElementCount = 1;
        sampleKeys = {"forward"};
        if (params[pi::PREPROCESS].as<bool>()) {
            dataFiles = filterDirContent(dataFiles, "_preproc.fastq");
        } else {
            dataFiles = filterDirContent(dataFiles, ".fastq");
        }
    } else if (subcall == pi::DETECT) {
        expectedElementCount = 1;
        sampleKeys = {"matched"};
        dataFiles = filterDirContent(dataFiles, "matched.sam");
    } else if (subcall == pi::ANALYZE) {
        expectedElementCount = 3;
        sampleKeys = {"splits", "singletonunassigned", "samplecounts"};
        dataFiles = filterDirContent(dataFiles, "_splits.sam", "singletonUnassigned.sam",
                                     "sampleCounts.tsv");
    }

    if (expectedElementCount != dataFiles.size()) {
        Logger::log(LogLevel::ERROR, "Expected ", expectedElementCount, " files, but found ",
                    dataFiles.size(), " files in ", conditionPath.string());
    }

    fs::path outConditionDir = fs::path(params["outdir"].as<std::string>()) /
                               params["subcall"].as<std::string>() / group /
                               conditionPath.filename();

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
// create ptree of the output files
pt::ptree Data::retrieveSampleOutputTree(fs::path outConditionDir, pt::ptree inputTree) {
    pt::ptree output;

    if (params["subcall"].as<std::string>() == pi::PREPROCESS) {
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
    } else if (params["subcall"].as<std::string>() == pi::ALIGN) {
        fs::path fwdInPath =
            fs::path(inputTree.get<std::string>("input.forward")).replace_extension(".sam");
        std::string fwdOutPath = replaceParentDirPath(outConditionDir, fwdInPath).string();
        std::string matched = addSuffix(fwdOutPath, "_matched", {});

        output.put("matched", matched);
    } else if (params["subcall"].as<std::string>() == pi::DETECT) {
        fs::path matchedInPath = fs::path(inputTree.get<std::string>("input.matched"));
        std::string splitsGeneralPath =
            replaceParentDirPath(outConditionDir, matchedInPath).string();
        std::string splitsOutPath = addSuffix(splitsGeneralPath, "splits", {"matched"});
        std::string multsplitsOutPath = addSuffix(splitsGeneralPath, "multsplits", {"matched"});
        std::string singletonUnassignedOutPath =
            addSuffix(splitsGeneralPath, "singletonUnassigned", {"matched"});
        std::string statsOutPath =
            fs::path(addSuffix(splitsGeneralPath, "sampleCounts", {"matched"}))
                .replace_extension(".tsv")
                .string();
        fs::path singletonTranscriptCountsOutPath =
            addSuffix(splitsGeneralPath, "singletonTranscriptCounts", {"matched"});
        singletonTranscriptCountsOutPath.replace_extension(".tsv");

        output.put("splits", splitsOutPath);
        output.put("multsplits", multsplitsOutPath);
        output.put("stats", statsOutPath);
        output.put("singleton", singletonTranscriptCountsOutPath.string());
        output.put("singletonUnassigned", singletonUnassignedOutPath);
    } else if (params["subcall"].as<std::string>() == pi::ANALYZE) {
        fs::path splits = fs::path(inputTree.get<std::string>("input.splits"));
        splits.replace_extension(".tab");
        std::string clu = replaceParentDirPath(outConditionDir, splits).string();
        std::string clusters = addSuffix(clu, "_clusters", {"_splits"});
        std::string clusterTranscriptCounts =
            addSuffix(clu, "_clusterTranscriptCounts", {"_splits"});
        fs::path supplementaryFeatures = addSuffix(clu, "_supplementaryFeatures", {"_splits"});
        supplementaryFeatures.replace_extension(".gff");
        output.put("clusters", clusters);
        output.put("clustertranscriptcounts", clusterTranscriptCounts);
        output.put("supplementaryfeatures", supplementaryFeatures.string());
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
    Logger::log(LogLevel::INFO, "Calling ", subcallStr, " pipeline");

    // create output directory (and subdirectory for subcall)
    createDirectories(outDir);         // create output
    createDirectories(outSubcallDir);  //  create subcall

    pt::ptree subcall;

    try {
        subcall = dataStructure.get_child(subcallStr);
    } catch (pt::ptree_error &e) {
        Logger::log(LogLevel::ERROR, subcallStr, " has not been found in the data structure");
    }

    std::deque<std::string> groups = {"trtms"};
    if (subcall.size() > 1) {
        groups.push_front("ctrls");
    }

    pt::ptree conditions, samples, path;
    fs::path outGroupDir, outConditionDir;

    for (unsigned i = 0; i < groups.size(); ++i) {
        // create directory for groups (e.g., ctrls, trtms)
        outGroupDir = outSubcallDir / fs::path(groups[i]);

        createDirectories(outGroupDir);

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

    Logger::log(LogLevel::INFO, "Finished ", subcallStr, " pipeline");
}

void Data::preproc() {
    // pipelines::preprocess::Preprocess srs(params);
    // callInAndOut(std::bind(&pipelines::preprocess::Preprocess::start, srs,
    // std::placeholders::_1));
}

void Data::align() {
    // Align aln(params);
    // callInAndOut(std::bind(&Align::start, aln, std::placeholders::_1));
}

void Data::splitReadCalling() {
    // Detect src(params);
    // callInAndOut(std::bind(&Detect::start, src, std::placeholders::_1));
}

void Data::analysis() {
    // Analyze clu(params);
    // callInAndOut(std::bind(&Analyze::start, &clu, std::placeholders::_1));
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

template <typename... Seeds>
PathVector Data::filterDirContent(const PathVector &paths, Seeds... seeds) {
    PathVector content;
    std::vector<std::string> seedsVec{seeds...};

    for (const auto &seed : seedsVec) {
        std::copy_if(paths.begin(), paths.end(), std::back_inserter(content), [&seed](fs::path x) {
            // only return if fs::path contains seed
            return x.string().find(seed) != std::string::npos;
        });
    }

    return content;
}

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
