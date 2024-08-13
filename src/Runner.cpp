#include "Runner.hpp"

#include <variant>

#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "DetectParameters.hpp"
#include "Logger.hpp"
#include "PreprocessParameters.hpp"

void Runner::runPipeline(int argc, const char *const *argv) {
    const auto parameters = ParameterParser::getParameters(argc, argv);

    std::visit(overloaded{[](auto) {
                              const auto message = "Invalid input";
                              Logger::log(LogLevel::ERROR, message);
                              throw std::runtime_error(message);
                          },
                          [](PreprocessParameters &params) { runPreprocessPipeline(params); },
                          [](AlignParameters &params) { runAlignPipeline(params); },
                          [](DetectParameters &params) { runDetectPipeline(params); },
                          [](AnalyzeParameters &params) { runAnalyzePipeline(params); },
                          [](CompleteParameters &params) { runCompletePipeline(params); }},
               parameters);
}
