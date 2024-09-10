#pragma once

// Boost
#include <boost/filesystem.hpp>

// Standard
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

// RNAnue
#include "DataTypes.hpp"
#include "Logger.hpp"

namespace annotation {

namespace fs = std::filesystem;

class FileType {
   public:
    enum Value : uint8_t { GFF, GTF };
    explicit constexpr FileType(Value p_value) : value(p_value) {};

    constexpr char attrDelim() const { return ';'; }

    constexpr char attributeAssignment() const {
        switch (value) {
            case GFF:
                return '=';
            case GTF:
                return ' ';
            default:
                return ' ';
        }
    }

    constexpr std::string defaultIDKey() const {
        switch (value) {
            case GFF:
                return "ID";
            case GTF:
                return "gene_id";
            default:
                return "gene_id";
        }
    }

    constexpr std::string defaultGroupKey() const {
        switch (value) {
            case GFF:
                return "Parent";
            case GTF:
                return "transcript_id";
            default:
                return "transcript_id";
        }
    }

    constexpr std::string defaultGeneNameKey() const { return "gene"; }

    constexpr operator Value() const { return value; }
    explicit operator bool() const = delete;

   private:
    Value value;
};

class FeatureParser {
   public:
    explicit FeatureParser(const std::unordered_set<std::string> &includedFeatures,
                           const std::optional<std::string> &featureIDFlag);
    ~FeatureParser() = default;

    dtp::FeatureMap parse(const fs::path featureFilePath) const;

   private:
    const std::unordered_set<std::string> includedFeatures;
    const std::optional<std::string> featureIDFlag;

    annotation::FileType getFileType(const fs::path &featureFilePath) const;
    dtp::FeatureMap iterateFeatureFile(const fs::path &featureFilePath,
                                       const annotation::FileType fileType) const;

    const std::unordered_map<std::string, std::string> getAttributes(
        const annotation::FileType fileType, const std::string &attributes) const;

    constexpr bool isValidFeature(const std::optional<std::vector<std::string>> &tokens) const;
};

}  // namespace annotation
