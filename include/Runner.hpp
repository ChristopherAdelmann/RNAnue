#pragma once

// Standard
#include <type_traits>
#include <variant>

// Classes
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "DetectParameters.hpp"
#include "Logger.hpp"
#include "ParameterParser.hpp"
#include "PreprocessParameters.hpp"

template <typename... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

class Runner {
   public:
    Runner() = delete;
    ~Runner() = delete;

    static void runPipeline(int argc, const char* const* argv);

   private:
    static void runPreprocessPipeline(const PreprocessParameters& preprocessParams);
    static void runAlignPipeline(const AlignParameters& alignParams);
    static void runDetectPipeline(const DetectParameters& detectParams);
    static void runAnalyzePipeline(const AnalyzeParameters& alignParams);

    static void runCompletePipeline(const CompleteParameters& completeParams);
};
