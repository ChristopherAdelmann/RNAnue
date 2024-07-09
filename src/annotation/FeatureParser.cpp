#include "FeatureParser.hpp"

using namespace Annotation;

/** @brief Constructs a FeatureParser object.
 *
 * This constructor initializes a FeatureParser object with the specified included features and
 * feature ID flag.
 *
 * @param includedFeatures The set of included features. Default is exon.
 *        Possible values include exon, gene, lncRNA, and other.
 * @param featureIDFlag The optional feature ID flag. Default is empty.
 *       Default is ID for GFF files and gene_id for GTF files.
 */
FeatureParser::FeatureParser(const std::unordered_set<std::string> &includedFeatures,
                             const std::optional<std::string> &featureIDFlag)
    : includedFeatures(includedFeatures), featureIDFlag(featureIDFlag) {}

dtp::FeatureMap FeatureParser::parse(const fs::path featureFilePath) const {
    FileType fileType = getFileType(featureFilePath);

    return iterateFeatureFile(featureFilePath, fileType);
}

FileType FeatureParser::getFileType(const fs::path &featureFilePath) const {
    std::ifstream file(featureFilePath.string());
    if (!file) {
        throw std::runtime_error("Could not open file: " + featureFilePath.string());
    }

    std::string line;
    std::getline(file, line);

    if (line.starts_with("##gff-version")) {
        return FileType(FileType::GFF);
    } else if (line.starts_with("##gtf-version")) {
        return FileType(FileType::GTF);
    }

    throw std::runtime_error(
        "Annotation file type not supported. First line does not contain "
        "##gff-version or ##gtf-version");
}

dtp::FeatureMap FeatureParser::iterateFeatureFile(const fs::path &featureFilePath,
                                                  const FileType fileType) const {
    dtp::FeatureMap featureMap;

    std::ifstream file(featureFilePath.string());

    if (!file.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", featureFilePath.string());
    }

    size_t parsedFeatures = 0;
    for (std::string line; std::getline(file, line);) {
        if (line[0] == '#') {
            continue;
        }

        auto getTokens =
            [](const std::string &line,
               const auto &includedFeatures) -> std::optional<std::vector<std::string>> {
            std::vector<std::string> tokens;
            std::istringstream issLine(line);
            for (std::string token; std::getline(issLine, token, '\t');) {
                // Checks at the third token if the feature is included
                if (tokens.size() == 2 && !includedFeatures.contains(token)) {
                    return std::nullopt;
                }

                tokens.push_back(token);
            }
            return tokens;
        };

        const auto tokens = getTokens(line, includedFeatures);

        if (!isValidFeature(tokens)) continue;

        const auto &tokens_v = tokens.value();
        const auto attributes = getAttributes(fileType, tokens_v[8]);

        auto getAttribute = [&attributes](const std::string &key) -> std::optional<std::string> {
            const auto it = attributes.find(key);
            if (it == attributes.end()) {
                return std::nullopt;
            }
            return it->second;
        };

        const std::string featureIDFlag = this->featureIDFlag.value_or(fileType.defaultIDKey());
        const auto identifier = getAttribute(featureIDFlag);

        if (!identifier.has_value()) {
            Logger::log(LogLevel::WARNING, "Could not find identifier in GFF file");
            continue;
        }

        const std::string &referenceID = tokens_v[0];
        const std::string &featureType = tokens_v[2];
        int startPosition = std::stoi(tokens_v[3]);
        int endPosition = std::stoi(tokens_v[4]);

        featureMap[referenceID].emplace_back(
        dtp::Feature{
            .referenceID = referenceID,
            .type = featureType,
            .startPosition = startPosition,
            .endPosition = endPosition,
            .strand = tokens_v[6][0] == '+' ? dtp::Strand::FORWARD : dtp::Strand::REVERSE,
            .id = identifier.value(),
            .groupID = getAttribute(fileType.defaultGroupKey())
            }
        );

        ++parsedFeatures;
    }

    const std::string includedFeatureTypes =
        std::accumulate(includedFeatures.begin(), includedFeatures.end(), std::string(),
                        [](const std::string &a, const std::string &b) { return a + b + ", "; });
    Logger::log(LogLevel::INFO, "Parsed ", std::to_string(parsedFeatures),
                                    " features of type: ", includedFeatureTypes);

    return featureMap;
}

std::unordered_map<std::string, std::string> const FeatureParser::getAttributes(
    const Annotation::FileType fileType, const std::string &attributes) const {
    std::unordered_map<std::string, std::string> attributeMap;
    std::istringstream issAttr(attributes);

    const char delim = fileType.attrDelim();

    for (std::string attribute; std::getline(issAttr, attribute, delim);) {
        const auto keyPosition = attribute.find(fileType.attributeAssignment());
        if (keyPosition == std::string::npos) continue;

        std::string key = attribute.substr(0, keyPosition);
        std::string value = attribute.substr(keyPosition + 1);

        if (key.empty() || value.empty()) continue;

        if (fileType == FileType::GTF) {
            const auto keyRet = std::ranges::remove(key, ' ');
            key.erase(keyRet.begin(), keyRet.end());
            const auto attributeRet = std::ranges::remove(value, '\"');
            attribute.erase(attributeRet.begin(), attributeRet.end());
        }

        attributeMap[key] = value;
    }

    return attributeMap;
}

constexpr bool FeatureParser::isValidFeature(
    const std::optional<std::vector<std::string>> &tokens) const {
    const std::array<char, 3> allowedStrand = {'+', '-', '.'};
    return tokens.has_value() && tokens.value().size() == 9 &&
           std::ranges::find(allowedStrand, tokens.value()[6][0]) != allowedStrand.end();
}