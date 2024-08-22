#include "Align.hpp"

#include <algorithm>
#include <variant>
#include <vector>

#include "AlignParameters.hpp"
#include "AlignSample.hpp"
#include "Constants.hpp"
#include "Logger.hpp"
#include "Utility.hpp"
#include "VariantOverload.hpp"

namespace pipelines {
namespace align {

void Align::process(const AlignData &data) {
    buildIndex();

    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    for (const auto &sample : data.treatmentSamples) {
        processSample(sample);
    }

    if (!data.controlSamples) {
        return;
    }

    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_CONTROL_MESSAGE);

    for (const auto &sample : *data.controlSamples) {
        processSample(sample);
    }
}

void Align::processSample(const AlignSampleType &sample) {
    buildIndex();

    std::visit(overloaded{[this](const AlignSampleSingle &sample) { processSingleEnd(sample); },
                          [this](const AlignSampleMergedPaired &sample) {
                              processMergedPairedEnd(sample);
                          }},
               sample);
}

void Align::processSingleEnd(const AlignSampleSingle &sample) {
    Logger::log(LogLevel::INFO, "Processing single end reads");

    alignSingleReads(sample.input.inputFastqPath, sample.output.outputAlignmentsPath);
    sortAlignmentsByQueryName(sample.output.outputAlignmentsPath,
                              sample.output.outputAlignmentsPath);
}

void Align::processMergedPairedEnd(const AlignSampleMergedPaired &sample) {
    Logger::log(LogLevel::INFO, "Processing merged paired end reads");

    alignSingleReads(sample.input.inputMergedFastqPath,
                     sample.output.outputAlignmentsMergedReadsPath);

    alignSingleReads(sample.input.inputSingletonForwardFastqPath,
                     sample.output.outputAlignmentsSingletonForwardReadsPath);

    alignSingleReads(sample.input.inputSingletonReverseFastqPath,
                     sample.output.outputAlignmentsSingletonReverseReadsPath);

    helper::mergeFiles(sample.output.outputAlignmentsPath,
                       {sample.output.outputAlignmentsMergedReadsPath,
                        sample.output.outputAlignmentsSingletonForwardReadsPath,
                        sample.output.outputAlignmentsSingletonReverseReadsPath});

    sortAlignmentsByQueryName(sample.output.outputAlignmentsPath,
                              sample.output.outputAlignmentsPath);
}

void Align::buildIndex() {
    fs::path referencePath = parameters.referenceGenome;
    int const threads = parameters.threadCount;

    fs::path indexFileName = referencePath.filename().replace_extension(".idx");
    indexPath = parameters.outputDir / indexFileName;

    if (fs::exists(indexPath)) {
        Logger::log(LogLevel::INFO, "Existing index found: ", indexPath);
        return;
    }

    Logger::log(LogLevel::INFO, "Building index");
    std::vector<std::string> args = {"-x", indexPath.string(),     "-d", referencePath.string(),
                                     "-t", std::to_string(threads)};

    auto c_args = convertToCStrings(args);

    int result = segemehl(c_args.size() - 1, c_args.data());

    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not create index for: ", referencePath);
    }
}

std::vector<std::string> Align::getGeneralAlignmentArgs() const {
    return {"-S",
            "-A",
            std::to_string(parameters.accuracy),
            "-U",
            std::to_string(parameters.minimumFragmentScore),
            "-W",
            std::to_string(parameters.minimumSpliceCoverage),
            "-Z",
            std::to_string(parameters.minimumFragmentLength),
            "-t",
            std::to_string(parameters.threadCount),
            "-m",
            std::to_string(parameters.minLengthThreshold),
            "-i",
            indexPath.string(),
            "-d",
            parameters.referenceGenome.string()};
}

void Align::alignSingleReads(const fs::path &queryFastqInPath,
                             const fs::path &alignmentsFastqOutPath) const {
    auto args = getGeneralAlignmentArgs();

    args.insert(args.end(),
                {"-q", queryFastqInPath.string(), "-o", alignmentsFastqOutPath.string()});

    auto c_args = convertToCStrings(args);

    int result = segemehl(c_args.size() - 1, c_args.data());

    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not align reads");
    }
}

void Align::alignPairedReads(const fs::path &queryForwardFastqInPath,
                             const fs::path &queryReverseFastqInPath,
                             const fs::path &alignmentsFastqOutPath) const {
    auto args = getGeneralAlignmentArgs();

    args.insert(args.end(),
                {"-q", queryForwardFastqInPath.string(), "-m", queryReverseFastqInPath.string(),
                 "-o", alignmentsFastqOutPath.string()});

    auto c_args = convertToCStrings(args);

    int result = segemehl(c_args.size() - 1, c_args.data());

    if (result != 0) {
        Logger::log(LogLevel::ERROR, "Could not align reads");
    }
}

void Align::sortAlignmentsByQueryName(const fs::path &alignmentsPath,
                                      const fs::path &sortedAlignmentsPath) const {
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
                                int(parameters.threadCount), &inFmt, &outFmt, emptyStr, true, 0);

    if (ret != 0) {
        Logger::log(LogLevel::ERROR, "Could not sort alignments");
    }

    Logger::log(LogLevel::INFO, "Sorting alignments done");
}

std::vector<char *> Align::convertToCStrings(std::vector<std::string> &args) const {
    std::vector<char *> c_args(args.size() + 1);
    std::transform(args.begin(), args.end(), c_args.begin(),
                   [](std::string &arg) { return const_cast<char *>(arg.c_str()); });
    c_args.back() = nullptr;  // argv must be null terminated

    return c_args;
}

}  // namespace align
}  // namespace pipelines
