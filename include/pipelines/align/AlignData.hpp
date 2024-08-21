#pragma once

// Standard
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

// Classes
#include "AlignSample.hpp"
#include "PreprocessData.hpp"
#include "pipelines/PipelineData.hpp"

namespace pipelines {
namespace align {

static const std::string validInputSuffixSingleton = preprocess::outSampleFastqSuffix;
static const std::array<std::string, 3> validInputSuffixesPaired{
    preprocess::outSampleFastqPairedMergeSuffix,
    preprocess::outSampleFastqPairedForwardSingletonSuffix,
    preprocess::outSampleFastqPairedReverseSingletonSuffix};

static const std::array<std::string, 4> validInputSuffixes{
    preprocess::outSampleFastqSuffix, preprocess::outSampleFastqPairedMergeSuffix,
    preprocess::outSampleFastqPairedForwardSingletonSuffix,
    preprocess::outSampleFastqPairedReverseSingletonSuffix};

static const std::string outSampleAlignedSuffix = "_aligned.bam";
static const std::string outSampleMergedAlignedSuffix = "_merged_aligned.bam";
static const std::string outSampleSingletonForwardAlignedSuffix = "_singleton_forward_aligned.bam";
static const std::string outSampleSingletonReverseAlignedSuffix = "_singleton_reverse_aligned.bam";

static constexpr std::string pipelinePrefix = "02_align";

struct AlignData : public pipelines::PipelineData {
    std::vector<AlignSampleType> treatmentSamples;
    std::optional<std::vector<AlignSampleType>> controlSamples;

    AlignData(const fs::path& outputDir, const fs::path& treatmentDir,
              const std::optional<fs::path> controlDir)
        : treatmentSamples(retrieveSamples(treatmentSampleGroup, treatmentDir, outputDir)),
          controlSamples(controlDir ? std::optional(retrieveSamples(controlSampleGroup,
                                                                    controlDir.value(), outputDir))
                                    : std::nullopt) {};

   private:
    static std::vector<AlignSampleType> retrieveSamples(const std::string& sampleGroup,
                                                        const fs::path& parentDir,
                                                        const fs::path& outputDir);
    static std::vector<InputSampleType> retrieveInputSamples(const fs::path& parentDir);
    static InputSampleType retrieveInputSample(const fs::path& sampleDir);
    static AlignInputPaired retrieveInputPaired(const std::string sampleName,
                                                const std::vector<fs::path>& inputSamples);
    static AlignInputSingle retrieveInputSingle(const std::string sampleName,
                                                const fs::path& inputSample);
};
}  // namespace align
}  // namespace pipelines
