#include "Runner.hpp"

#include <filesystem>
#include <optional>
#include <variant>

#include "Align.hpp"
#include "AlignData.hpp"
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "DetectParameters.hpp"
#include "Logger.hpp"
#include "Preprocess.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "pipelines/PipelineData.hpp"

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

    const auto pipeline = preprocess::Preprocess(parameters);
    pipeline.process(data);
}

void Runner::runAlignPipeline(const AlignParameters &parameters) {
    Logger::log(LogLevel::INFO, "Running align pipeline");

    const auto inputDirs = InputDirectories(parameters.outputDir, preprocess::pipelinePrefix);

    const auto data = align::AlignData(parameters.outputDir, inputDirs.treatmentInputDir,
                                       inputDirs.controlInputDir);

    auto pipeline = align::Align(parameters);
    pipeline.process(data);
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

    runPreprocessPipeline(parameters.preprocessParameters);
    runAlignPipeline(parameters.alignParameters);
    runDetectPipeline(parameters.detectParameters);
    runAnalyzePipeline(parameters.analyzeParameters);

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
