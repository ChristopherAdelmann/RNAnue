#include "Utility.hpp"

#include <filesystem>
#include <fstream>

#include "Logger.hpp"
#include "seqan3/io/sam_file/input.hpp"

namespace helper {
void createTmpDir(fs::path path) {
    Logger::log(LogLevel::INFO, "Create temporary directory " + path.string());
    deleteDir(path);
    fs::create_directory(path);
}

// delete folder
void deleteDir(fs::path path) {
    if (fs::exists(path)) {
        fs::remove_all(path);
    }
}

std::optional<fs::path> getDirIfExists(const fs::path& path) {
    if (fs::is_directory(path)) {
        return path;
    }

    return std::nullopt;
}

std::string getUUID() {
    boost::uuids::random_generator uuidGenerator;
    return boost::uuids::to_string(uuidGenerator());
}

void mergeSamFiles(std::vector<fs::path> inputPaths, fs::path outputPath) {
    if (inputPaths.empty()) {
        Logger::log(LogLevel::WARNING, "No input files to merge");
        return;
    }

    seqan3::sam_file_output outputFile{outputPath};

    for (const auto& inputPath : inputPaths) {
        seqan3::sam_file_input inputFile{inputPath};
        inputFile | outputFile;
    }
}

void concatAndDeleteFilesInTmpDir(const fs::path& tmpDir, const fs::path& outPath) {
    std::ofstream outStream(outPath, std::ios::binary);

    for (const auto& entry : fs::directory_iterator(tmpDir)) {
        std::ifstream inStream(entry.path(), std::ios::binary);
        outStream << inStream.rdbuf();
        inStream.close();
    }
}

bool isContained(const int32_t value, const int32_t comparisonValue, const int32_t tolerance) {
    return value >= comparisonValue - tolerance && value <= comparisonValue + tolerance;
}

std::size_t countUniqueSamEntries(fs::path path) {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::unordered_set<std::string> uniqueIDs;
    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty() && line[0] != '@') {
            std::string id = line.substr(0, line.find('\t'));
            uniqueIDs.insert(id);
        }
    }

    return uniqueIDs.size();
}

std::size_t countSamEntriesSeqAn(fs::path path) {
    seqan3::sam_file_input fin{path.string(), seqan3::fields<>{}};
    return std::ranges::distance(fin.begin(), fin.end());
}

size_t countSamEntries(fs::path path) {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    int count = 0;
    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty() && line[0] != '@') {
            ++count;
        }
    }

    return count;
}

std::string generateRandomHexColor() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, 255);

    std::stringstream ss;
    ss << "#";
    for (int i = 0; i < 3; ++i) {
        ss << std::setfill('0') << std::setw(2) << std::hex << distr(gen);
    }

    return ss.str();  // Return the hex color code as a string
}

std::string getTime() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t current_time = std::chrono::system_clock::to_time_t(now);

    std::ostringstream time_stream;
    time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]") << " ";

    return time_stream.str();
}

Timer::Timer() : start(std::chrono::high_resolution_clock::now()) {}
Timer::~Timer() { stop(); }
void Timer::stop() {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}

}  // namespace helper
