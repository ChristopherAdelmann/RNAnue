#include "Utility.hpp"

void helper::createTmpDir(fs::path path) {
    Logger::log(LogLevel::INFO, "Create temporary directory " + path.string());
    deleteDir(path);
    fs::create_directory(path);
}

// delete folder
void helper::deleteDir(fs::path path) {
    if (fs::exists(path)) {
        fs::remove_all(path);
    }
}

void helper::mergeSamFiles(std::vector<fs::path> inputPaths, fs::path outputPath) {
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

std::size_t helper::countUniqueSamEntries(fs::path path) {
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

std::size_t helper::countSamEntriesSeqAn(fs::path path) {
    seqan3::sam_file_input fin{path.string(), seqan3::fields<>{}};
    return std::ranges::distance(fin.begin(), fin.end());
}

size_t helper::countSamEntries(fs::path path) {
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

std::string helper::generateRandomHexColor() {
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

void helper::printTree(const boost::property_tree::ptree& pt, int level = 0) {
    for (const auto& node : pt) {
        std::cout << std::string(level * 2, ' ') << node.first << ": "
                  << node.second.get_value<std::string>() << "\n";
        printTree(node.second, level + 1);
    }
}

std::string helper::getTime() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t current_time = std::chrono::system_clock::to_time_t(now);

    std::ostringstream time_stream;
    time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]") << " ";

    return time_stream.str();
}

helper::Timer::Timer() : start(std::chrono::high_resolution_clock::now()) {}
helper::Timer::~Timer() { stop(); }
void helper::Timer::stop() {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}