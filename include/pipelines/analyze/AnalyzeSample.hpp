#pragma once

// Standard
#include <filesystem>
#include <string>

namespace pipelines {
namespace analyze {
namespace fs = std::filesystem;

struct AnalyzeInput {
    std::string sampleName;

    fs::path splitAlignmentsPath;
    fs::path multisplitAlignmentsPath;
    fs::path unassignedContiguousAlignmentsPath;

    fs::path contiguousAlignmentsTranscriptCountsPath;

    fs::path sampleFragmentCountsPath;
};

struct AnalyzeOutput {
    fs::path interactionsPath;
    fs::path interactionsTranscriptCountsPath;

    fs::path interactionsBEDPath;
    fs::path interactionsBEDARCPath;

    fs::path supplementaryFeaturesPath;
};

struct AnalyzeSample {
    AnalyzeInput input;
    AnalyzeOutput output;
};

}  // namespace analyze
}  // namespace pipelines
