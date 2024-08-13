#pragma once

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "Logger.hpp"
#include "PreprocessSample.hpp"
#include "pipelines/PipelineData.hpp"

namespace pipelines {
namespace preprocess {

namespace fs = std::filesystem;

static constexpr std::array<std::string, 4> validSuffixes{".fastq", ".fastq.gz", ".fq", ".fq.gz"};
static constexpr std::array<std::string, 3> validForwardTags{"_R1", "_1", "_forward"};
static constexpr std::array<std::string, 3> validReverseTags{"_R2", "_2", "_reverse"};
static constexpr std::string pipelinePrefix = "01_preprocess";

struct PreprocessData : public pipelines::PipelineData {
    std::vector<PreprocessSampleType> treatmentSamples;
    std::optional<std::vector<PreprocessSampleType>> controlSamples;

    PreprocessData(const fs::path& outputDir, const fs::path& treatmentDir,
                   const std::optional<fs::path> controlDir)
        : treatmentSamples(retrieveSamples(treatmentDir, outputDir)),
          controlSamples(controlDir ? std::optional(retrieveSamples(controlDir.value(), outputDir))
                                    : std::nullopt) {};

   private:
    static std::vector<PreprocessSampleType> retrieveSamples(const fs::path& parentDir,
                                                             const fs::path& outputDir);
    static std::vector<InputSampleType> retrieveInputSamples(const fs::path& parentDir);
    static InputSampleType retrieveInputSample(const fs::path& sampleDir);
    static bool validatePairedFilePaths(std::pair<fs::path, fs::path>& pairedInputPaths);
};
}  // namespace preprocess
}  // namespace pipelines
