#include "DetectData.hpp"

namespace pipelines {
namespace detect {

const std::vector<DetectSample> DetectData::retrieveSamples(const std::string& sampleGroup,
                                                            const fs::path& parentDir,
                                                            const fs::path& outputDir) {
    const std::vector<InputSample> inputSamples = retrieveInputSamples(parentDir);

    std::vector<DetectSample> samples;
    samples.reserve(inputSamples.size());

    const fs::path outputDirPipeline = outputDir / sampleGroup / pipelinePrefix;

    for (const InputSample& inputSample : inputSamples) {
        const fs::path outputDirSample = outputDirPipeline / inputSample.sampleName;

        const fs::path outputSplitAlignmentsPath =
            outputDirSample / (inputSample.sampleName + outSampleSplitAlignmentsSuffix);
        const fs::path outputMultisplitAlignmentsPath =
            outputDirSample / (inputSample.sampleName + outSampleMultisplitAlignmentsSuffix);
        const fs::path outputUnassignedContiguousAlignmentsPath =
            outputDirSample /
            (inputSample.sampleName + outSampleUnassignedContiguousAlignmentsSuffix);
        const fs::path outputContiguousAlignmentsTranscriptCountsPath =
            outputDirSample /
            (inputSample.sampleName + outSampleContiguousAlignmentsTranscriptCountsSuffix);

        const fs::path outputSharedStatsPath = outputDirPipeline / outSharedStatsSuffix;

        samples.push_back(DetectSample(
            inputSample,
            OutputSample{outputSplitAlignmentsPath, outputMultisplitAlignmentsPath,
                         outputUnassignedContiguousAlignmentsPath,
                         outputContiguousAlignmentsTranscriptCountsPath, outputSharedStatsPath}));
    }

    return samples;
}

const std::vector<InputSample> DetectData::retrieveInputSamples(const fs::path& parentDir) {
    const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

    std::vector<InputSample> samples;

    for (const fs::path& sampleDir : sampleDirs) {
        const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, {validInputSuffix});

        if (sampleFiles.size() != 1) {
            const std::string message = "Expected 1 input file in " + sampleDir.string() +
                                        " but found " + std::to_string(sampleFiles.size());
            Logger::log(LogLevel::ERROR, message);
            throw std::runtime_error(message);
        }

        const fs::path& sampleFile = sampleFiles.front();
        const auto sampleName = getSampleName(sampleFile);

        if (!hasSuffix(sampleFile, validInputSuffix)) {
            const std::string message = "Invalid input file suffix: " + sampleFile.string() +
                                        " . Expected: " + validInputSuffix;
            Logger::log(LogLevel::ERROR, message);
            throw std::runtime_error(message);
        }

        samples.emplace_back(sampleName, sampleFile);
    }

    return samples;
}

}  // namespace detect
}  // namespace pipelines
