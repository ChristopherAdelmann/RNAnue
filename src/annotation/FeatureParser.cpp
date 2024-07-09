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
        throw std::runtime_error("Could not open file: " + featureFilePath.string());
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

        const std::string featureIDFlag = this->featureIDFlag.value_or(fileType.featureIDFlag());
        const auto identifier = getIdentifier(fileType, tokens_v[8], featureIDFlag);

        if (!identifier.has_value()) {
            Logger::log(LogLevel::WARNING, "Could not find identifier in GFF file");
            continue;
        }

        const std::string seqid = tokens.value()[0];

        dtp::Feature feature{
            .referenceID = seqid,
            .type = tokens_v[2],
            .startPosition = std::stoi(tokens_v[3]),
            .endPosition = std::stoi(tokens_v[4]),
            .strand = tokens_v[6][0] == '+' ? dtp::Strand::FORWARD : dtp::Strand::REVERSE,
            .id = identifier.value()};

        featureMap[seqid].push_back(feature);

        ++parsedFeatures;
    }

    const std::string includedFeatureTypes =
        std::accumulate(includedFeatures.begin(), includedFeatures.end(), std::string(),
                        [](const std::string &a, const std::string &b) { return a + b + ", "; });
    Logger::log(LogLevel::INFO, "Parsed " + std::to_string(parsedFeatures) +
                                    " features of type:" + includedFeatureTypes);

    return featureMap;
}

std::optional<std::string> FeatureParser::getIdentifier(const FileType fileType,
                                                        const std::string &attributes,
                                                        const std::string &query) const {
    const char delim = fileType.attrDelim();
    std::istringstream issAttr(attributes);
    for (std::string attr; std::getline(issAttr, attr, delim);) {
        const size_t matchPos = attr.find(query);
        if (matchPos == std::string::npos) continue;

        std::string id = attr.substr(matchPos + query.size() + 1);

        switch (fileType) {
            case FileType::GFF:
                return id;
            case FileType::GTF:
                const auto ret = std::ranges::remove(id, '\"');
                id.erase(ret.begin(), ret.end());
                return id;
        }
    }
    return std::nullopt;
}

constexpr bool FeatureParser::isValidFeature(
    const std::optional<std::vector<std::string>> &tokens) const {
    const std::array<char, 3> allowedStrand = {'+', '-', '.'};
    return tokens.has_value() && tokens.value().size() == 9 &&
           std::ranges::find(allowedStrand, tokens.value()[6][0]) != allowedStrand.end();
}