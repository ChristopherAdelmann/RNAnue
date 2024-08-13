#pragma once

// Standard
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <variant>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Classes
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "Closing.hpp"
#include "CompleteParameters.hpp"
#include "Config.hpp"
#include "Constants.hpp"
#include "DetectParameters.hpp"
#include "GeneralParameters.hpp"
#include "Logger.hpp"
#include "ParameterOptions.hpp"
#include "PreprocessParameters.hpp"

namespace po = boost::program_options;

class ParameterParser {
   public:
    using ParametersVariant = std::variant<CompleteParameters, PreprocessParameters,
                                           AlignParameters, DetectParameters, AnalyzeParameters>;

    static ParametersVariant getParameters(int argc, const char* const argv[]);

    ParameterParser() = delete;

   private:
    static po::variables_map parseParameters(int argc, const char* const argv[]);

    static void insertConfigFileParameters(po::variables_map& params);

    static po::options_description getCommandLineOptions();
    static po::options_description getConfigFileOptions();

    static po::positional_options_description getPositionalOptions();

    static void printVersion();
};
