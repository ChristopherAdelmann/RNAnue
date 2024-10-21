#pragma once

// Standard
#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>

// RNAnue
#include "FileType.hpp"
#include "GenomicFeature.hpp"

namespace annotation {

namespace fs = std::filesystem;

class FeatureParser {
   public:
    explicit FeatureParser(const std::unordered_set<std::string> &includedFeatures,
                           const std::optional<std::string> &featureIDFlag);
    FeatureParser(const FeatureParser &) = default;
    FeatureParser(FeatureParser &&) = delete;
    auto operator=(const FeatureParser &) -> FeatureParser & = delete;
    auto operator=(FeatureParser &&) -> FeatureParser & = delete;
    ~FeatureParser() = default;

    [[nodiscard]] auto parse(const fs::path &featureFilePath) const -> dataTypes::FeatureMap;

   private:
    std::unordered_set<std::string> includedFeatures;
    std::optional<std::string> featureIDFlag;

    [[nodiscard]] auto iterateFeatureFile(const fs::path &featureFilePath,
                                          annotation::FileType fileType) const
        -> dataTypes::FeatureMap;

    [[nodiscard]] static auto getFileType(const fs::path &featureFilePath) -> annotation::FileType;
    [[nodiscard]] static auto getTokens(const std::string &line,
                                        const std::unordered_set<std::string> &includedFeatures)
        -> std::optional<std::vector<std::string>>;
    [[nodiscard]] static auto getAttributes(annotation::FileType fileType,
                                            const std::string &attributes)
        -> std::unordered_map<std::string, std::string>;

    [[nodiscard]] static constexpr auto isValidFeature(
        const std::optional<std::vector<std::string>> &tokens) -> bool;
};

}  // namespace annotation
