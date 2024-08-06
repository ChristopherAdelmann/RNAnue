#pragma once

// Boost
#include <boost/filesystem.hpp>

// Standard
#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

// RNAnue
#include "DataTypes.hpp"
#include "Logger.hpp"

namespace Annotation {
namespace fs = boost::filesystem;

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

    Annotation::FileType getFileType(const fs::path &featureFilePath) const;
    dtp::FeatureMap iterateFeatureFile(const fs::path &featureFilePath,
                                       const Annotation::FileType fileType) const;

    const std::unordered_map<std::string, std::string> getAttributes(
        const Annotation::FileType fileType, const std::string &attributes) const;

    constexpr bool isValidFeature(const std::optional<std::vector<std::string>> &tokens) const;
};

}  // namespace Annotation
