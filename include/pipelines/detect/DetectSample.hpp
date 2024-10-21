#pragma once

// Standard
#include <filesystem>
#include <string>

namespace pipelines::detect {
namespace fs = std::filesystem;

struct DetectInput {
    std::string sampleName;
    fs::path inputAlignmentsPath;
};

struct DetectOutput {
    fs::path outputSplitAlignmentsPath;
    fs::path outputMultisplitAlignmentsPath;
    fs::path outputUnassignedContiguousAlignmentsPath;
    fs::path outputContiguousAlignmentsTranscriptCountsPath;

    fs::path outputSharedReadCountsPath;
};

struct DetectSample {
    DetectInput input;
    DetectOutput output;
};

}  // namespace pipelines::detect
