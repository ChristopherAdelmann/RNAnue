#pragma once

// Standard
#include <filesystem>
#include <variant>

namespace pipelines {
namespace preprocess {

namespace fs = std::filesystem;

struct PreprocessInput {
    std::string sampleName;
};

struct PreprocessSampleInputSingle : public PreprocessInput {
    fs::path inputFastqPath;

    PreprocessSampleInputSingle(const std::string& sampleName, const fs::path& samplePath)
        : PreprocessInput{sampleName}, inputFastqPath{samplePath} {};
};

struct PreprocessSampleInputPaired : public PreprocessInput {
    fs::path inputForwardFastqPath;
    fs::path inputReverseFastqPath;

    PreprocessSampleInputPaired(const std::string& sampleName, const fs::path& forwardPath,
                                const fs::path& reversePath)
        : PreprocessInput{sampleName},
          inputForwardFastqPath{forwardPath},
          inputReverseFastqPath{reversePath} {};
};

using InputSampleType = std::variant<PreprocessSampleInputSingle, PreprocessSampleInputPaired>;

struct PrepocessSampleOutputSingle {
    fs::path outputFastqPath;
};

struct PrepocessSampleOutputPaired {
    fs::path outputMergedFastqPath;
    fs::path outputSingletonForwardFastqPath;
    fs::path outputSingletonReverseFastqPath;
};

using OutputSampleType = std::variant<PrepocessSampleOutputSingle, PrepocessSampleOutputPaired>;

struct PreprocessSampleSingle {
    PreprocessSampleInputSingle input;
    PrepocessSampleOutputSingle output;
};

struct PreprocessSamplePaired {
    PreprocessSampleInputPaired input;
    PrepocessSampleOutputPaired output;
};

using PreprocessSampleType = std::variant<PreprocessSampleSingle, PreprocessSamplePaired>;

}  // namespace preprocess
}  // namespace pipelines
