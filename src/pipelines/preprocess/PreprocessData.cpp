#include "pipelines/preprocess/PreprocessData.hpp"

#include "Logger.hpp"
#include "PreprocessSample.hpp"

namespace pipelines {
namespace preprocess {
std::vector<PreprocessSampleType> PreprocessData::retrieveSamples(const fs::path& parentDir,
                                                                  const fs::path& outputDir) {
    const std::vector<InputSampleType> inputSamples = retrieveInputSamples(parentDir);

    std::vector<PreprocessSampleType> samples;
    samples.reserve(inputSamples.size());

    for (const InputSampleType& inputSample : inputSamples) {
        if (const auto* inputSampleSingle =
                std::get_if<PreprocessSampleInputSingle>(&inputSample)) {
            const std::string parentName = inputSampleSingle->inputFastqPath.parent_path().stem();
            const fs::path outputSampleFastqPath =
                outputDir / parentName / (inputSampleSingle->sampleName + "_passed.fastq.gz");

            samples.push_back(PreprocessSampleSingle{*inputSampleSingle, {outputSampleFastqPath}});

            continue;

        } else if (const auto* inputSamplePaired =
                       std::get_if<PreprocessSampleInputPaired>(&inputSample)) {
            const std::string parentName =
                inputSamplePaired->inputForwardFastqPath.parent_path().stem();
            const fs::path outputDirSample = outputDir / parentName;

            const fs::path outputSampleFastqPathMerged =
                outputDirSample / (inputSamplePaired->sampleName + "_merged_passed.fastq.gz");
            const fs::path outputSampleFastqPathForwardSingleton =
                outputDirSample /
                (inputSamplePaired->sampleName + "_forward_singleton_passed_R2.fastq.gz");
            const fs::path outputSampleFastqPathReverseSingleton =
                outputDirSample /
                (inputSamplePaired->sampleName + "_reverse_singleton_passed_R1.fastq.gz");

            samples.push_back(PreprocessSamplePaired{
                *inputSamplePaired,
                {outputSampleFastqPathMerged, outputSampleFastqPathForwardSingleton,
                 outputSampleFastqPathReverseSingleton}});
            continue;
        }
    }

    return samples;
}

std::vector<InputSampleType> PreprocessData::retrieveInputSamples(const fs::path& parentDir) {
    const std::vector<fs::path> sampleDirs = getDirectories(parentDir);

    std::vector<InputSampleType> samples;

    for (const fs::path& sampleDir : sampleDirs) {
        samples.push_back(retrieveInputSample(sampleDir));
    }

    return samples;
}

InputSampleType PreprocessData::retrieveInputSample(const fs::path& sampleDir) {
    const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, validSuffixes);
    const auto numSamples = sampleFiles.size();

    using namespace std::string_literals;
    if (numSamples == 0 || numSamples > 2) {
        const std::string message = (numSamples == 0 ? "No"s : "More than two"s) +
                                    " valid sample files found in " + sampleDir.string();
        Logger::log(LogLevel::ERROR, message);
        throw std::runtime_error(message);
    }

    const auto& firstFile = sampleFiles.front();
    const auto sampleName = getSampleName(firstFile);

    if (numSamples == 1) {
        if (hasAnySuffix(firstFile.stem(), validForwardTags) ||
            hasAnySuffix(firstFile.stem(), validReverseTags)) {
            const auto message =
                "Found single valid file, but expected two based on file ending in " +
                sampleDir.string();
            Logger::log(LogLevel::ERROR, message);
            throw std::runtime_error(message);
        }
        return PreprocessSampleInputSingle{sampleName, firstFile};
    }

    std::pair<fs::path, fs::path> pairedInputPaths{firstFile, sampleFiles.back()};

    if (!validatePairedFilePaths(pairedInputPaths)) {
        const auto message = "Invalid paired file paths found in " + sampleDir.string();
        Logger::log(LogLevel::ERROR, message);
        throw std::runtime_error(message);
    }

    return PreprocessSampleInputPaired{sampleName, pairedInputPaths.first, pairedInputPaths.second};
}

bool PreprocessData::validatePairedFilePaths(std::pair<fs::path, fs::path>& pairedInputPaths) {
    const bool validForward = hasAnySuffix(pairedInputPaths.first.stem(), validForwardTags);

    if (!validForward) {
        pairedInputPaths = {pairedInputPaths.second, pairedInputPaths.first};

        const bool validForward = hasAnySuffix(pairedInputPaths.first.stem(), validForwardTags);

        if (!validForward) {
            return false;
        }
    }

    const bool validReverse = hasAnySuffix(pairedInputPaths.second.stem(), validReverseTags);

    return validReverse;
}
}  // namespace preprocess
}  // namespace pipelines
