#include "Utility.hpp"

#include <filesystem>
#include <fstream>

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

void mergeFiles(const fs::path& outputPath, const std::vector<fs::path>& paths) {
    std::ofstream outfile(outputPath, std::ios::binary);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    for (const auto& path : paths) {
        std::ifstream inputFilestream(path, std::ios::binary);

        if (!inputFilestream.is_open()) {
            throw std::runtime_error("Could not open file");
        }

        outfile << inputFilestream.rdbuf();
    }
}

void mergeSamFiles(std::vector<fs::path> inputPaths, fs::path outputPath) {
    if (inputPaths.empty()) {
        return;
    }

    std::ofstream outfile(outputPath);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    // Only add the header from the first file
    std::ifstream firstFile(inputPaths[0]);
    std::string line;
    while (std::getline(firstFile, line)) {
        if (line.empty() || line[0] == '@') {
            outfile << line << "\n";
        } else {
            break;
        }
    }

    // Merge the contents of all the files
    for (const auto& path : inputPaths) {
        std::ifstream inputFile(path);
        while (std::getline(inputFile, line)) {
            if (line.empty() || line[0] == '@') {
                continue;  // Skip the header lines
            }
            outfile << line << "\n";
        }
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

void printTree(const boost::property_tree::ptree& pt, int level = 0) {
    for (const auto& node : pt) {
        std::cout << std::string(level * 2, ' ') << node.first << ": "
                  << node.second.get_value<std::string>() << "\n";
        printTree(node.second, level + 1);
    }
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
