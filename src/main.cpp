// Boost
#include <boost/program_options.hpp>

// Standard
#include <execinfo.h>
#include <unistd.h>

#include <array>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

// Internal

#include "Config.hpp"
#include "ParameterParser.hpp"
#include "Runner.hpp"

void showVersion() {
    const std::string versionString =
        "RNAnue v" + std::to_string(RNAnue_VERSION_MAJOR) + "." +
        std::to_string(RNAnue_VERSION_MINOR) + "." + std::to_string(RNAnue_VERSION_PATCH) + " - " +
        "Detect RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.";

    Logger::log(LogLevel::INFO, versionString);
}

// A helper function to simplify the main part.
template <class T>
auto operator<<(std::ostream& ostream, const std::vector<T>& vec) -> std::ostream& {
    copy(vec.begin(), vec.end(), std::ostream_iterator<T>(ostream, " "));
    return ostream;
}

namespace std {
auto operator<<(std::ostream& ostream, const std::vector<std::string>& vec) -> std::ostream& {
    for (const auto& item : vec) {
        ostream << item << " ";
    }
    return ostream;
}
}  // namespace std

void handler(int sig) {
    constexpr size_t MAX_FRAMES = 10;
    std::array<void*, MAX_FRAMES> array{};
    int size = backtrace(array.data(), MAX_FRAMES);

    // print out all the frames to stderr
    std::cerr << "Error: signal " << sig << ":" << std::endl;
    backtrace_symbols_fd(array.data(), size, STDERR_FILENO);
    exit(1);
}

auto main(int argc, const char* const argv[]) -> int {
    signal(SIGSEGV, handler);

    Runner::runPipeline(argc, argv);
}
