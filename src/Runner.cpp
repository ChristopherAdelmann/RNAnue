#include "Runner.hpp"

#include <variant>

#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "DetectParameters.hpp"
#include "Logger.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"

void Runner::runPipeline(int argc, const char *const argv[]) {
    const auto parameters = ParameterParser::getParameters(argc, argv);

    std::visit(Pipeline(), parameters);
}

void Runner::runPreprocessPipeline(const preprocess::PreprocessParameters &parameters) {
    Logger::log(LogLevel::INFO, "Running preprocess pipeline");

    if (!parameters.preprocessEnabled) {
        Logger::log(LogLevel::INFO, "Preprocess pipeline is disabled in the parameters");
        return;
    }

    const auto data = preprocess::PreprocessData(parameters.outputDir, parameters.treatmentsDir,
                                                 parameters.controlDir);
}

void Runner::runAlignPipeline(const AlignParameters &parameters) {
    Logger::log(LogLevel::INFO, "Running align pipeline");

    Logger::log(LogLevel::ERROR, "Align pipeline is not implemented yet");
}

void Runner::runDetectPipeline(const DetectParameters &parameters) {
    Logger::log(LogLevel::INFO, "Running detect pipeline");

    Logger::log(LogLevel::ERROR, "Detect pipeline is not implemented yet");
}

void Runner::runAnalyzePipeline(const AnalyzeParameters &parameters) {
    Logger::log(LogLevel::INFO, "Running analyze pipeline");

    Logger::log(LogLevel::ERROR, "Analyze pipeline is not implemented yet");
}

void Runner::runCompletePipeline(const CompleteParameters &parameters) {
    Logger::log(LogLevel::INFO, "Running complete pipeline");

    Logger::log(LogLevel::ERROR, "Complete pipeline is not implemented yet");
}
// TODO: Implement the rest of the pipeline functions

void Runner::Pipeline::operator()(const preprocess::PreprocessParameters &params) {
    Runner::runPreprocessPipeline(params);
};
void Runner::Pipeline::operator()(const AlignParameters &params) {
    Runner::runAlignPipeline(params);
};
void Runner::Pipeline::operator()(const DetectParameters &params) {
    Runner::runDetectPipeline(params);
};
void Runner::Pipeline::operator()(const AnalyzeParameters &params) {
    Runner::runAnalyzePipeline(params);
};
void Runner::Pipeline::operator()(const CompleteParameters params) {
    Runner::runCompletePipeline(params);
};
