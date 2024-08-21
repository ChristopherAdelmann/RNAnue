#pragma once

// Standard
#include <filesystem>
#include <string>
#include <variant>

namespace pipelines {
namespace align {
namespace fs = std::filesystem;

struct AlignInput {
    std::string sampleName;
};

struct AlignInputSingle : public AlignInput {
    fs::path inputFastq;

    AlignInputSingle(const std::string& sampleName, const fs::path& inputFastq)
        : AlignInput{sampleName}, inputFastq{inputFastq} {}
};

struct AlignInputPaired : public AlignInput {
    fs::path inputMergedFastqPath;
    fs::path inputSingletonForwardFastqPath;
    fs::path inputSingletonReverseFastqPath;

    AlignInputPaired(const std::string& sampleName, const fs::path& inputMergedFastqPath,
                     const fs::path& inputSingletonForwardFastqPath,
                     const fs::path& inputSingletonReverseFastqPath)
        : AlignInput{sampleName},
          inputMergedFastqPath{inputMergedFastqPath},
          inputSingletonForwardFastqPath{inputSingletonForwardFastqPath},
          inputSingletonReverseFastqPath{inputSingletonReverseFastqPath} {}
};

using InputSampleType = std::variant<AlignInputSingle, AlignInputPaired>;

struct AlignOutputSingle {
    fs::path outputAlignmentsPath;
};

struct AlignOutputPaired {
    fs::path outputAlignmentsPath;
    fs::path outputAlignmentsMergedReadsPath;
    fs::path outputAlignmentsSingletonForwardReadsPath;
    fs::path outputAlignmentsSingletonReverseReadsPath;
};

using OutputSampleType = std::variant<AlignOutputSingle, AlignOutputPaired>;

struct AlignSampleSingle {
    AlignInputSingle input;
    AlignOutputSingle output;
};

struct AlignSamplePaired {
    AlignInputPaired input;
    AlignOutputPaired output;
};

using AlignSampleType = std::variant<AlignSampleSingle, AlignSamplePaired>;

}  // namespace align
}  // namespace pipelines
