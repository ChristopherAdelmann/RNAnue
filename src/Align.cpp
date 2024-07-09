#include "Align.hpp"

Align::Align(po::variables_map params) : params(params) {}

void Align::alignReads(const std::string &query, const std::string &mate,
                       const std::string &matched) {
    Logger::log(LogLevel::INFO, "Aligning reads");

    std::vector<std::string> args = {"-S",
                                     "-A",
                                     std::to_string(params["accuracy"].as<int>()),
                                     "-U",
                                     std::to_string(params["minfragsco"].as<int>()),
                                     "-W",
                                     std::to_string(params["minsplicecov"].as<int>()),
                                     "-Z",
                                     std::to_string(params["minfraglen"].as<int>()),
                                     "-t",
                                     std::to_string(params["threads"].as<int>()),
                                     "-m",
                                     std::to_string(params["minlen"].as<std::size_t>()),
                                     "-i",
                                     indexPath.string(),
                                     "-d",
                                     params["dbref"].as<std::string>(),
                                     "-q",
                                     query};

    if (!mate.empty()) {
        args.insert(args.end(), {"-p", mate});
    }

    args.insert(args.end(), {"-o", matched});

    std::vector<char *> c_args(args.size() + 1);
    std::transform(args.begin(), args.end(), c_args.begin(),
                   [](const std::string &arg) { return const_cast<char *>(arg.c_str()); });
    c_args.back() = nullptr;  // argv must be null terminated

    int result = segemehl(c_args.size() - 1, c_args.data());

    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not align reads");
    }
}

void Align::buildIndex() {
    fs::path referencePath = params["dbref"].as<std::string>();
    int const threads = params["threads"].as<int>();

    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    fs::path indexFileName = referencePath.filename().replace_extension(".idx");
    indexPath = outDir / indexFileName;

    if (fs::exists(indexPath)) {
        Logger::log(LogLevel::INFO, "Existing index found: ", indexPath);
        return;
    }

    Logger::log(LogLevel::INFO, "Building index");
    std::vector<std::string> args = {"-x", indexPath.string(),     "-d", referencePath.string(),
                                     "-t", std::to_string(threads)};

    std::vector<char *> c_args(args.size() + 1);
    std::transform(args.begin(), args.end(), c_args.begin(),
                   [](std::string &arg) { return const_cast<char *>(arg.c_str()); });
    c_args.back() = nullptr;  // argv must be null terminated

    int result = segemehl(c_args.size() - 1, c_args.data());

    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not create index for: ", referencePath);
    }
}

void Align::sortAlignmentsByQueryName(const std::string &alignmentsPath,
                                      const std::string &sortedAlignmentsPath) {
    Logger::log(LogLevel::INFO, "Sorting alignments");

    // TODO Adapt output format to selection from config
    const size_t SORT_DEFAULT_MEGS_PER_THREAD = 768;
    const size_t maxMem = SORT_DEFAULT_MEGS_PER_THREAD << 20;
    const htsFormat inFmt = {sequence_data, sam, {1, 6}, no_compression, 0, 0};
    const htsFormat outFmt = {sequence_data, sam, {1, 6}, no_compression, 0, 0};

    const fs::path tempDir = fs::path(alignmentsPath).parent_path();
    char emptyStr[] = "";
    const char wbStr[] = "wb";

    int ret = bam_sort_core_ext(QueryName, emptyStr, 0, true, true, alignmentsPath.c_str(),
                                tempDir.c_str(), sortedAlignmentsPath.c_str(), wbStr, maxMem,
                                params["threads"].as<int>(), &inFmt, &outFmt, emptyStr, true, 0);

    if (ret != 0) {
        Logger::log(LogLevel::ERROR, "Could not sort alignments");
    }

    Logger::log(LogLevel::INFO, "Sorting alignments done");
}

void Align::start(pt::ptree sample) {
    buildIndex();

    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    // align the reads (regular)
    std::string inFwd = input.get<std::string>("forward");
    std::string outMatched = output.get<std::string>("matched");

    alignReads(inFwd, "", outMatched);

    sortAlignmentsByQueryName(outMatched, outMatched);

    // TODO implement unpaired read analysis with all downstream steps
    // if (params["readtype"].as<std::string>() == "PE") {
    //     if (params["unprd"].as<bool>()) {
    //         std::string inR1only = input.get<std::string>("R1only");
    //         std::string outR1only = output.get<std::string>("matched_R1only");
    //         alignReads(inR1only, "", outR1only);

    //         std::string inR2only = input.get<std::string>("R2only");
    //         std::string outR2only = output.get<std::string>("matched_R2only");
    //         alignReads(inR2only, "", outR2only);
    //     }

    //     if (params["unmrg"].as<bool>()) {
    //         std::string inR1unmrg = input.get<std::string>("R1unmerged");
    //         std::string inR2unmrg = input.get<std::string>("R2unmerged");
    //         std::string outunmrg = output.get<std::string>("matched_unmerged");
    //         alignReads(inR1unmrg, inR2unmrg, outunmrg);
    //     }
    // }
}