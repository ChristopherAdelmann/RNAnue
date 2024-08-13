#include "ParameterParser.hpp"

#include "ParameterOptions.hpp"

ParameterParser::ParametersVariant ParameterParser::getParameters(int argc,
                                                                  const char *const *argv) {
    const auto params = parseParameters(argc, argv);

    const std::string subcall = params["subcall"].as<std::string>();

    if (subcall == constants::pipelines::COMPLETE) {
        return CompleteParameters{params};
    } else if (subcall == constants::pipelines::PREPROCESS) {
        return PreprocessParameters{params};
    } else if (subcall == constants::pipelines::ALIGN) {
        return AlignParameters{params};
    } else if (subcall == constants::pipelines::DETECT) {
        return DetectParameters{params};
    } else if (subcall == constants::pipelines::ANALYZE) {
        return AnalyzeParameters{params};
    } else {
        Logger::log(LogLevel::ERROR, "Unknown subcall: " + subcall);
        exit(EXIT_FAILURE);
    }
}

po::variables_map ParameterParser::parseParameters(int argc, const char *const argv[]) {
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

    if (params.count("version")) {
        Closing::printQuote();
        exit(EXIT_SUCCESS);
    }

    if (params.count("help")) {
        std::cout << commandLineOptions << std::endl;
        Closing::printQuote();
        exit(EXIT_SUCCESS);
    }

    if (!params.count("subcall")) {
        Logger::log(LogLevel::ERROR, "Please provide a subcall.");
    }

    Logger::setLogLevel(params["loglevel"].as<std::string>());

    insertConfigFileParameters(params);

    return params;
}

void ParameterParser::insertConfigFileParameters(po::variables_map &params) {
    if (params.count("config") == 0) return;

    const po::options_description configFileOptions{getConfigFileOptions()};

    const std::string configFilePath{params["config"].as<std::string>()};

    std::ifstream configIn{configFilePath};

    if (!configIn) {
        Logger::log(LogLevel::ERROR, "Configuration file could not be opened!");
    }

    po::store(po::parse_config_file(configIn, configFileOptions), params);
    notify(params);
}

po::options_description ParameterParser::getCommandLineOptions() {
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

po::options_description ParameterParser::getConfigFileOptions() {
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

po::positional_options_description ParameterParser::getPositionalOptions() {
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
