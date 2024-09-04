#pragma once

// Standard
#include <filesystem>
#include <optional>
#include <type_traits>
#include <variant>

// Classes
#include "Align.hpp"
#include "AlignData.hpp"
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "Detect.hpp"
#include "DetectData.hpp"
#include "DetectParameters.hpp"
#include "Logger.hpp"
#include "ParameterParser.hpp"
#include "Preprocess.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "Utility.hpp"
#include "pipelines/PipelineData.hpp"

using namespace pipelines;

class Runner {
   public:
    Runner() = delete;
    ~Runner() = delete;

    static void runPipeline(int argc, const char* const argv[]);

   private:
    struct InputDirectories {
        InputDirectories(const fs::path parentDir, const std::string pipelinePrefix)
            : treatmentInputDir(parentDir / pipelinePrefix / treatmentSampleGroup),
              controlInputDir(
                  helper::getDirIfExists(parentDir / pipelinePrefix / controlSampleGroup)) {};

        fs::path treatmentInputDir;
        std::optional<fs::path> controlInputDir;
    };

    struct Pipeline {
        void operator()(const preprocess::PreprocessParameters& params);
        void operator()(const AlignParameters& params);
        void operator()(const DetectParameters& params);
        void operator()(const AnalyzeParameters& params);
        void operator()(const CompleteParameters params);
    };

    static void runPreprocessPipeline(const preprocess::PreprocessParameters& preprocessParams);
    static void runAlignPipeline(const AlignParameters& alignParams);
    static void runDetectPipeline(const DetectParameters& detectParams);
    static void runAnalyzePipeline(const AnalyzeParameters& alignParams);

    static void runCompletePipeline(const CompleteParameters& completeParams);
};
