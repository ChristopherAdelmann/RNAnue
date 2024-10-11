#pragma once

// Standard
#include <array>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

// Classes
#include "Logger.hpp"
#include "Utility.hpp"

namespace pipelines {

namespace fs = std::filesystem;

using namespace helper;

template <typename Container>
concept StringContainer = std::ranges::range<Container> &&
                          std::same_as<std::ranges::range_value_t<Container>, std::string>;

static const std::string treatmentSampleGroup = "treatment";
static const std::string controlSampleGroup = "control";

struct PipelineData {
   protected:
    static bool isHidden(const std::filesystem::path &path) {
        std::string filename = path.filename().string();
        return !filename.empty() && filename[0] == '.';
    }

    static std::string getSampleName(const fs::path &filePath) {
        std::string fileName = filePath.stem().string();
        size_t underscorePos = fileName.find("_");
        if (underscorePos != std::string::npos) {
            return fileName.substr(0, underscorePos);
        } else {
            return fileName;
        }
    }

    static bool contains(const std::string &fullString, const std::string &substring) {
        return fullString.find(substring) != std::string::npos;
    };

    template <StringContainer Container>
    static bool containsAny(const std::string &fullString, const Container &substrings) {
        for (const auto &substring : substrings) {
            if (contains(fullString, substring)) {
                return true;
            }
        }
        return false;
    };

    static std::vector<fs::path> getSubDirectories(const fs::path &parentDir) {
        std::vector<fs::path> directories;
        for (const auto &entry : fs::directory_iterator(parentDir)) {
            if (entry.is_directory()) {
                directories.push_back(entry.path());
            } else if (entry.is_regular_file() && !isHidden(entry)) {
                Logger::log(LogLevel::WARNING,
                            "Found file in directory, which is not allowed: ", entry.path());
            }
        }
        return directories;
    };
};
}  // namespace pipelines
