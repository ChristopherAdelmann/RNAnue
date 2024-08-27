#pragma once

// Standard
#include <filesystem>
#include <string>
#include <vector>

namespace pipelines {
namespace detect {
namespace fs = std::filesystem;

struct InputSample {
    std::string sampleName;
    fs::path inputAlignmentsPath;
};

struct OutputSample {
    fs::path outputSplitAlignmentsPath;
    fs::path outputMultisplitAlignmentsPath;
    fs::path outputUnassignedContiguousAlignmentsPath;
    fs::path outputContiguousAlignmentsTranscriptCountsPath;

    fs::path outputSharedStatsPath;
};

struct DetectSample {
    InputSample input;
    OutputSample output;
};

}  // namespace detect
}  // namespace pipelines
