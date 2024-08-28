#include "AlignData.hpp"

#include <variant>
#include <vector>

#include "AlignSample.hpp"

namespace pipelines {
namespace align {
std::vector<AlignSampleType> AlignData::retrieveSamples(const std::string& sampleGroup,
                                                        const fs::path& parentDir,
                                                        const fs::path& outputDir) {
    const std::vector<InputSampleType> inputSamples = retrieveInputSamples(parentDir);

    std::vector<AlignSampleType> samples;
    samples.reserve(inputSamples.size());

    const fs::path outputDirPipeline = outputDir / pipelinePrefix / sampleGroup;

    for (const InputSampleType& inputSample : inputSamples) {
        if (const auto* inputSampleSingle = std::get_if<AlignInputSingle>(&inputSample)) {
            const fs::path outputDirSample = outputDirPipeline / inputSampleSingle->sampleName;

            fs::create_directories(outputDirSample);

            const fs::path outputAlignmentsPath =
                outputDirPipeline / inputSampleSingle->sampleName /
                (inputSampleSingle->sampleName + outSampleAlignedSuffix);

            samples.push_back(AlignSampleSingle{*inputSampleSingle, {outputAlignmentsPath}});

            const auto message = "Single-end sample " + inputSampleSingle->sampleName + " found";
            Logger::log(LogLevel::INFO, message);

            continue;

        } else if (const auto* inputSamplePaired = std::get_if<AlignInputPaired>(&inputSample)) {
            const std::string parentName = inputSamplePaired->sampleName;
            const fs::path outputDirSample = outputDirPipeline / inputSamplePaired->sampleName;

            fs::create_directories(outputDirSample);

            const fs::path outputAlignmentsPath =
                outputDirSample / (parentName + outSampleAlignedSuffix);

            const fs::path outputAlignmentsMergedReadsPath =
                outputDirSample / (parentName + outSampleMergedAlignedSuffix);
            const fs::path outputAlignmentsSingletonForwardPath =
                outputDirSample / (parentName + outSampleSingletonForwardAlignedSuffix);
            const fs::path outputAlignmentsSingletonReversePath =
                outputDirSample / (parentName + outSampleSingletonReverseAlignedSuffix);

            samples.push_back(AlignSampleMergedPaired{
                *inputSamplePaired,
                {outputAlignmentsPath, outputAlignmentsMergedReadsPath,
                 outputAlignmentsSingletonForwardPath, outputAlignmentsSingletonReversePath}});

            const auto message = "Paired-end sample " + parentName + " found";
            Logger::log(LogLevel::INFO, message);

            continue;
        }
    }

    return samples;
}

std::vector<InputSampleType> AlignData::retrieveInputSamples(const fs::path& parentDir) {
    const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

    std::vector<InputSampleType> samples;

    for (const fs::path& sampleDir : sampleDirs) {
        samples.push_back(retrieveInputSample(sampleDir));
    }

    return samples;
}

InputSampleType AlignData::retrieveInputSample(const fs::path& sampleDir) {
    const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, validInputSuffixes);
    const auto numSamples = sampleFiles.size();

    using namespace std::string_literals;
    if (numSamples != validInputSuffixSingleton.size() ||
        numSamples != validInputSuffixesPaired.size()) {
        const std::string message =
            "Found invalid number of sample files (" + std::to_string(numSamples) + ")" +
            ". Expectected either " + std::to_string(validInputSuffixSingleton.size()) + " or " +
            std::to_string(validInputSuffixesPaired.size()) + " in " + sampleDir.string();
        Logger::log(LogLevel::ERROR, message);
        throw std::runtime_error(message);
    }

    const auto& firstFile = sampleFiles.front();
    const auto sampleName = getSampleName(firstFile);

    // First check paired because singleton suffix is subset of paired suffix
    if (numSamples == validInputSuffixesPaired.size()) {
        return retrieveInputPaired(sampleName, sampleFiles);
    }

    return retrieveInputSingle(sampleName, firstFile);
}

AlignInputPaired AlignData::retrieveInputPaired(const std::string sampleName,
                                                const std::vector<fs::path>& inputSamples) {
    std::optional<fs::path> mergedFastqPath;
    std::optional<fs::path> singletonForwardFastqPath;
    std::optional<fs::path> singletonReverseFastqPath;

    for (const fs::path& path : inputSamples) {
        if (hasSuffix(path, preprocess::outSampleFastqPairedMergeSuffix)) {
            mergedFastqPath = path;
        } else if (hasSuffix(path, preprocess::outSampleFastqPairedForwardSingletonSuffix)) {
            singletonForwardFastqPath = path;
        } else if (hasSuffix(path, preprocess::outSampleFastqPairedReverseSingletonSuffix)) {
            singletonReverseFastqPath = path;
        } else {
            const std::string message =
                "Unexpected file " + path.string() + " found in " + path.parent_path().string();
            Logger::log(LogLevel::WARNING, message);
        }
    }

    if (!mergedFastqPath.has_value() || !singletonForwardFastqPath.has_value() ||
        !singletonReverseFastqPath.has_value()) {
        const std::string message = "The directory " + inputSamples.front().parent_path().string() +
                                    " is missing one or more of the following files: " +
                                    preprocess::outSampleFastqPairedMergeSuffix + ", " +
                                    preprocess::outSampleFastqPairedForwardSingletonSuffix + ", " +
                                    preprocess::outSampleFastqPairedReverseSingletonSuffix;
        Logger::log(LogLevel::ERROR, message);
        throw std::runtime_error(message);
    }

    return AlignInputPaired{sampleName, mergedFastqPath.value(), singletonForwardFastqPath.value(),
                            singletonReverseFastqPath.value()};
}

AlignInputSingle AlignData::retrieveInputSingle(const std::string sampleName,
                                                const fs::path& inputSample) {
    if (!hasSuffix(inputSample, preprocess::outSampleFastqSuffix)) {
        const std::string message = "The directory " + inputSample.parent_path().string() +
                                    " is missing the file: " + preprocess::outSampleFastqSuffix;
        Logger::log(LogLevel::ERROR, message);
        throw std::runtime_error(message);
    }
    return AlignInputSingle{sampleName, inputSample};
}

}  // namespace align
}  // namespace pipelines
