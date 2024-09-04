#pragma once

// Standard
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

// Classes
#include "AnalyzeSample.hpp"
#include "DetectData.hpp"
#include "pipelines/PipelineData.hpp"

namespace pipelines {
namespace analyze {
namespace fs = std::filesystem;

static const std::string validInputSplitAlignmentsSuffix = detect::outSampleSplitAlignmentsSuffix;
static const std::string validInputMultisplitAlignmentsSuffix =
    detect::outSampleMultisplitAlignmentsSuffix;
static const std::string validInputUnassignedContiguousAlignmentsSuffix =
    detect::outSampleUnassignedContiguousAlignmentsSuffix;
static const std::string validInputContiguousAlignmentsTranscriptCountsSuffix =
    detect::outSampleContiguousAlignmentsTranscriptCountsSuffix;
static const std::string validSharedReadCountsSuffix = detect::outSampleReadCountsSummarySuffix;

static const std::array<std::string, 5> validSuffices = {
    validInputSplitAlignmentsSuffix, validInputMultisplitAlignmentsSuffix,
    validInputUnassignedContiguousAlignmentsSuffix,
    validInputContiguousAlignmentsTranscriptCountsSuffix, validSharedReadCountsSuffix};

static const std::string outInteractionsSuffix = "_interactions.tsv";
static const std::string outInteractionsTranscrptCountsSuffix =
    "_interaction_transcript_counts.tsv";
static const std::string outInteractionsBEDSuffix = "_interaction_regions.bed";
static const std::string outInteractionsBEDARCSuffix = "_interaction_regions.arc";
static const std::string outSupplementaryFeaturesSuffix = "_supplementary_features.gff";

static constexpr std::string pipelinePrefix = "04_analyze";

struct AnalyzeData : public PipelineData {
    std::vector<AnalyzeSample> treatmentSamples;
    std::optional<std::vector<AnalyzeSample>> controlSamples;

    AnalyzeData(const fs::path& outputDir, const fs::path& treatmentDir,
                const std::optional<fs::path> controlDir)
        : treatmentSamples(retrieveSamples(treatmentSampleGroup, treatmentDir, outputDir)),
          controlSamples(controlDir ? std::optional<std::vector<AnalyzeSample>>(retrieveSamples(
                                          controlSampleGroup, controlDir.value(), outputDir))
                                    : std::nullopt) {}

   private:
    static const std::vector<AnalyzeSample> retrieveSamples(const std::string& sampleGroup,
                                                            const fs::path& parentDir,
                                                            const fs::path& outputDir);

    static const std::vector<AnalyzeInput> retrieveInputSamples(const fs::path& parentDir);
};

}  // namespace analyze
}  // namespace pipelines
