#include "ParameterParser.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "Closing.hpp"
#include "Constants.hpp"
#include "Logger.hpp"
#include "ParameterOptions.hpp"

namespace pipelines {
auto ParameterParser::getParameters(int argc, const char *const argv[])  // NOLINT
    -> ParameterParser::ParametersVariant {
    const auto params = parseParameters(argc, argv);

    const std::string subcall = params["subcall"].as<std::string>();

    if (subcall == constants::pipelines::COMPLETE) {
        return CompleteParameters{params};
    }
    if (subcall == constants::pipelines::PREPROCESS) {
        return preprocess::PreprocessParameters{params};
    }
    if (subcall == constants::pipelines::ALIGN) {
        return AlignParameters{params};
    }
    if (subcall == constants::pipelines::DETECT) {
        return DetectParameters{params};
    }
    if (subcall == constants::pipelines::ANALYZE) {
        return AnalyzeParameters{params};
    }

    Logger::log(LogLevel::ERROR, "Unknown subcall: " + subcall);
    exit(EXIT_FAILURE);
}

auto ParameterParser::parseParameters(int argc,
                                      const char *const argv[]) -> po::variables_map {  // NOLINT
    const po::options_description commandLineOptions{getCommandLineOptions()};

    const po::positional_options_description positionalOptions{getPositionalOptions()};

    po::variables_map params;
    store(po::command_line_parser(argc, argv)
              .options(commandLineOptions)
              .positional(positionalOptions)
              .run(),
          params);

    notify(params);

    printVersion();

    if (params.count("version") != 0U) {
        Closing::printQuote();
        exit(EXIT_SUCCESS);
    }

    if (params.count("help") != 0U) {
        std::cout << commandLineOptions << std::endl;
        Closing::printQuote();
        exit(EXIT_SUCCESS);
    }

    if (params.count("subcall") == 0U) {
        Logger::log(LogLevel::ERROR, "Please provide a subcall.");
    }

    Logger::setLogLevel(params["loglevel"].as<std::string>());

    insertConfigFileParameters(params);

    return params;
}

void ParameterParser::insertConfigFileParameters(po::variables_map &params) {
    if (params.count("config") == 0) {
        return;
    }

    const po::options_description configFileOptions{getConfigFileOptions()};

    const std::string configFilePath{params["config"].as<std::string>()};

    std::ifstream configIn{configFilePath};

    if (!configIn) {
        Logger::log(LogLevel::ERROR, "Configuration file could not be opened!");
    }

    po::store(po::parse_config_file(configIn, configFileOptions), params);
    notify(params);
}

auto ParameterParser::getCommandLineOptions() -> po::options_description {
    const po::options_description generalOptions{ParameterOptions::getGeneralOptions()};
    const po::options_description preprocessOptions{ParameterOptions::getPreprocessOptions()};
    const po::options_description alignOptions{ParameterOptions::getAlignOptions()};
    const po::options_description detectOptions{ParameterOptions::getDetectOptions()};
    const po::options_description analyzeOptions{ParameterOptions::getAnalyzeOptions()};
    const po::options_description otherOptions{ParameterOptions::getOtherOptions()};
    const po::options_description subcallOptions{ParameterOptions::getSubcallOptions()};

    po::options_description commandLineOptions{"Command line options"};

    commandLineOptions.add(generalOptions)
        .add(preprocessOptions)
        .add(alignOptions)
        .add(detectOptions)
        .add(analyzeOptions)
        .add(otherOptions)
        .add(subcallOptions);

    return commandLineOptions;
}

auto ParameterParser::getConfigFileOptions() -> po::options_description {
    const po::options_description generalOptions{ParameterOptions::getGeneralOptions()};
    const po::options_description preprocessOptions{ParameterOptions::getPreprocessOptions()};
    const po::options_description alignOptions{ParameterOptions::getAlignOptions()};
    const po::options_description detectOptions{ParameterOptions::getDetectOptions()};
    const po::options_description analyzeOptions{ParameterOptions::getAnalyzeOptions()};

    po::options_description configFileOptions{"Config file options"};

    configFileOptions.add(generalOptions)
        .add(preprocessOptions)
        .add(alignOptions)
        .add(detectOptions)
        .add(analyzeOptions);

    return configFileOptions;
}

auto ParameterParser::getPositionalOptions() -> po::positional_options_description {
    po::positional_options_description positionalOptions;
    positionalOptions.add("subcall", 1);

    return positionalOptions;
}

void ParameterParser::printVersion() {
    const std::string versionString =
        "RNAnue v" + std::to_string(RNAnue_VERSION_MAJOR) + "." +
        std::to_string(RNAnue_VERSION_MINOR) + "." + std::to_string(RNAnue_VERSION_PATCH) + " - " +
        "Detect RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.";

    Logger::log(LogLevel::INFO, versionString);
}

}  // namespace pipelines
