// Boost
#include <boost/program_options.hpp>

// Standard
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <array>
#include <iostream>
#include <string>

// Classes

#include "Base.hpp"
#include "Closing.hpp"
#include "Config.hpp"
#include "Constants.hpp"
#include "FeatureAnnotator.hpp"
#include "ParameterParser.hpp"
#include "Runner.hpp"

namespace po = boost::program_options;

void showVersion() {
    const std::string versionString =
        "RNAnue v" + std::to_string(RNAnue_VERSION_MAJOR) + "." +
        std::to_string(RNAnue_VERSION_MINOR) + "." + std::to_string(RNAnue_VERSION_PATCH) + " - " +
        "Detect RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.";

    Logger::log(LogLevel::INFO, versionString);
}

// A helper function to simplify the main part.
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

namespace std {
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vec) {
    for (const auto& item : vec) {
        os << item << " ";
    }
    return os;
}
}  // namespace std

void handler(int sig) {
    std::array<void*, 10> array;
    size_t size = backtrace(array.data(), 10);

    // print out all the frames to stderr
    std::cerr << "Error: signal " << sig << ":" << std::endl;
    backtrace_symbols_fd(array.data(), size, STDERR_FILENO);
    exit(1);
}

int main(int argc, const char* const argv[]) {
    signal(SIGSEGV, handler);

    Runner::runPipeline(argc, argv);
}
