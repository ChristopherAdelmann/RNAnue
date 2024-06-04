#include "Align.hpp"

Align::Align(po::variables_map params)
    : params(params), segemehlSysCall(params["segemehl"].as<std::string>()) {}

void Align::alignReads(const std::string &query, const std::string &mate,
                       const std::string &matched) {
    Logger::log(LogLevel::INFO, "Aligning reads");
    std::string align = segemehlSysCall;
    align += " -S ";  // split mode
    align += " -A " + std::to_string(params["accuracy"].as<int>());
    align += " -U " + std::to_string(params["minfragsco"].as<int>());
    align += " -W " + std::to_string(params["minsplicecov"].as<int>());
    align += " -Z " + std::to_string(params["minfraglen"].as<int>());
    align += " -t " + std::to_string(params["threads"].as<int>());
    align += " -m " + std::to_string(params["minlen"].as<std::size_t>());
    align += " -i " + indexPath.string();
    align += " -d " + params["dbref"].as<std::string>();
    align += " -q " + query;
    if (!mate.empty()) {
        align += " -p " + mate;
    }
    align += " -o " + matched;
    const char *alignCallChar = align.c_str();
    int result = system(alignCallChar);
    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not align reads");
        exit(1);
    }
}

void Align::buildIndex() {
    fs::path referencePath = params["dbref"].as<std::string>();
    int const threads = params["threads"].as<int>();

    fs::path outDir = fs::path(params["outdir"].as<std::string>());
    indexPath = outDir / referencePath.replace_extension(".idx").filename();

    if (fs::exists(indexPath)) {
        Logger::log(LogLevel::INFO, "Existing index found: ", indexPath);
        return;
    }

    Logger::log(LogLevel::INFO, "Building index");
    std::string generateIndexCall = segemehlSysCall + " -x " + indexPath.string() + "-t" +
                                    std::to_string(threads) + " -d " + referencePath.string();
    int result = system(generateIndexCall.c_str());
    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not create index for: ", referencePath);
        exit(1);
    }
}

void Align::sortAlignments(const std::string &alignmentsPath) {
    Logger::log(LogLevel::INFO, "Sorting alignments");

    const std::string samtoolsSortCall = "samtools sort -n -@ " +
                                         std::to_string(params["threads"].as<int>()) + " -o " +
                                         alignmentsPath + " " + alignmentsPath;

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

    // sortAlignments(outMatched);

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