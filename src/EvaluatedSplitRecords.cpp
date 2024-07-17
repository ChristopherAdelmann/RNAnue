#include "EvaluatedSplitRecords.hpp"

// TODO Needs optimization and implementation for multiple splits
std::optional<EvaluatedSplitRecords> EvaluatedSplitRecords::calculateEvaluatedSplitRecords(
    SplitRecords &splitRecords, const BaseParameters &parameters) {
    if (splitRecords.size() != 2) {
        Logger::log(LogLevel::DEBUG, "Currently only two split records are supported!");
        return std::nullopt;
    }

    const seqan3::dna5_vector &seq1 = splitRecords[0].sequence();
    const seqan3::dna5_vector &seq2 = splitRecords[1].sequence();

    std::optional<CoOptimalPairwiseAligner::AlignmentResult> complementarity =
        calculateComplementarity(seq1, seq2, parameters.minComplementarity,
                                 parameters.minComplementarityFraction);

    if (!complementarity.has_value()) {
        return std::nullopt;
    }

    std::optional<HybridizationResult> hybridization =
        calculateHybridization(seq1, seq2, parameters.mfeThreshold);

    if (!hybridization.has_value()) {
        return std::nullopt;
    }

    addTagsToRecords(splitRecords, complementarity.value(), hybridization.value());

    return EvaluatedSplitRecords{splitRecords, complementarity.value(), hybridization.value()};
}

std::optional<EvaluatedSplitRecords> EvaluatedSplitRecords::calculateEvaluatedSplitRecords(
    SplitRecords &splitRecords, const SplicingParameters &parameters,
    const Annotation::FeatureAnnotator &featureAnnotator) {
    if (splitRecords.size() != 2) {
        Logger::log(LogLevel::DEBUG, "Currently only two split records are supported!");
        return std::nullopt;
    }

    if (!splitRecords[0].reference_position().has_value() ||
        !splitRecords[1].reference_position().has_value()) [[unlikely]] {
        Logger::log(LogLevel::WARNING, "Could not determine reference position of split records: ",
                    splitRecords[0].id(), " or ", splitRecords[1].id());
        return std::nullopt;
    }

    if (splitRecords[0].reference_id() != splitRecords[1].reference_id()) {
        return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
    }

    const auto &record1 =
        splitRecords[0].reference_position().value() < splitRecords[1].reference_position().value()
            ? splitRecords[0]
            : splitRecords[1];
    const auto &record2 =
        splitRecords[0].reference_position().value() < splitRecords[1].reference_position().value()
            ? splitRecords[1]
            : splitRecords[0];

    const auto featureRecord1 = featureAnnotator.getBestOverlappingFeature(
        record1, parameters.referenceIDs, parameters.orientation);

    if (!featureRecord1.has_value() || !featureRecord1.value().groupID.has_value()) {
        return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
    }

    const auto featureRecord2 = featureAnnotator.getBestOverlappingFeature(
        record2, parameters.referenceIDs, parameters.orientation);

    if (!featureRecord2.has_value() || !featureRecord2.value().groupID.has_value()) {
        return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
    }

    if (featureRecord1.value().groupID != featureRecord2.value().groupID) {
        return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
    }

    // Check if record1 is at splice junction start
    const auto record1EndPosition = dtp::recordEndPosition(record1).value();

    bool record1AtSpliceJunctionStart = helper::isContained(
        record1EndPosition, featureRecord1.value().endPosition - 1, parameters.splicingTolerance);

    if (!record1AtSpliceJunctionStart) {
        return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
    }

    const auto record2StartPosition = record2.reference_position().value();
    bool record2AtSpliceJunctionEnd = helper::isContained(
        record2StartPosition, featureRecord2.value().startPosition, parameters.splicingTolerance);

    if (!record2AtSpliceJunctionEnd) {
        return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
    }

    const std::string &referenceID = parameters.referenceIDs[record1.reference_id().value()];

    const auto inBetweenRegion =
        dtp::GenomicRegion{referenceID, record1EndPosition, record2StartPosition};

    auto it = featureAnnotator.overlappingFeatureIterator(inBetweenRegion, parameters.orientation);

    const bool hasInBetweenExon =
        std::any_of(it.begin(), it.end(), [&featureRecord1](const auto &feature) {
            return feature.groupID.has_value() &&
                   feature.groupID.value() == featureRecord1.value().groupID;
        });

    if (!hasInBetweenExon) {
        return std::nullopt;
    }

    return calculateEvaluatedSplitRecords(splitRecords, parameters.baseParameters);
}

// TODO Implement overall score calculation instead of comparing individual scores
bool EvaluatedSplitRecords::operator<(const EvaluatedSplitRecords &other) const {
    // Calculate if fraction difference is more or equal than 0.05
    if (std::abs(complementarityResult.fraction - other.complementarityResult.fraction) >= 0.05) {
        return complementarityResult.fraction < other.complementarityResult.fraction;
    }

    // Calculate if complementarity difference is more or equal than 0.05
    if (std::abs(complementarityResult.complementarity -
                 other.complementarityResult.complementarity) >= 0.05) {
        return complementarityResult.complementarity < other.complementarityResult.complementarity;
    }

    // Calculate if energy difference is more or equal than 0.05
    if (std::abs(hybridizationResult.energy - other.hybridizationResult.energy) >= 0.05) {
        return hybridizationResult.energy < other.hybridizationResult.energy;
    }

    // If both have crosslinking results prefer the one with higher crosslinking score
    if (hybridizationResult.crosslinkingResult.has_value() &&
        other.hybridizationResult.crosslinkingResult.has_value()) {
        return hybridizationResult.crosslinkingResult.value().normCrosslinkingScore <
               other.hybridizationResult.crosslinkingResult.value().normCrosslinkingScore;
    }

    return false;
}

bool EvaluatedSplitRecords::operator>(const EvaluatedSplitRecords &other) const {
    return !(*this < other);
}

void EvaluatedSplitRecords::addTagsToRecords(
    SplitRecords &splitRecords, const CoOptimalPairwiseAligner::AlignmentResult &complementarity,
    const HybridizationResult &hybridization) {
    for (auto &record : splitRecords) {
        // Complementarity tags
        const int length =
            complementarity.endPositions.first - complementarity.beginPositions.first;
        record.tags()["XL"_tag] = length;
        record.tags()["XC"_tag] = static_cast<float>(complementarity.complementarity);
        record.tags()["XR"_tag] = static_cast<float>(complementarity.fraction);
        record.tags()["XS"_tag] = complementarity.score;
        // rec1.tags()["XA"_tag] = res.a;
        // rec2.tags()["XA"_tag] = res.b;
        // rec1.tags()["XM"_tag] = res.matches;
        // rec2.tags()["XM"_tag] = res.matches;

        // Hybridization tags
        record.tags()["XE"_tag] = static_cast<float>(hybridization.energy);
        if (hybridization.crosslinkingResult.has_value()) {
            record.tags()["XN"_tag] =
                static_cast<float>(hybridization.crosslinkingResult.value().normCrosslinkingScore);
            record.tags()["XP"_tag] =
                hybridization.crosslinkingResult.value().preferredCrosslinkingScore;
            record.tags()["XQ"_tag] =
                hybridization.crosslinkingResult.value().nonPreferredCrosslinkingScore;
            record.tags()["XW"_tag] =
                hybridization.crosslinkingResult.value().wobbleCrosslinkingScore;
        }
    }
}

std::optional<CoOptimalPairwiseAligner::AlignmentResult>
EvaluatedSplitRecords::calculateComplementarity(const seqan3::dna5_vector &seq1,
                                                const seqan3::dna5_vector &seq2,
                                                const double minComplementarity,
                                                const double minComplementarityFraction) {
    const auto seq1Reverse = seq1 | std::views::reverse;

    const auto results = CoOptimalPairwiseAligner{complementaryScoringScheme()}.getLocalAlignments(
        std::tie(seq1Reverse, seq2));
    return getOptimalAlignment(results, minComplementarity, minComplementarityFraction);
}

constexpr seqan3::nucleotide_scoring_scheme<int8_t>
EvaluatedSplitRecords::complementaryScoringScheme() {
    seqan3::nucleotide_scoring_scheme scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}};

    scheme.score('A'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'A'_dna5) = 1;
    scheme.score('G'_dna5, 'C'_dna5) = 1;
    scheme.score('C'_dna5, 'G'_dna5) = 1;
    scheme.score('G'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'G'_dna5) = 1;

    scheme.score('T'_dna5, 'T'_dna5) = -1;
    scheme.score('A'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'G'_dna5) = -1;

    scheme.score('A'_dna5, 'G'_dna5) = -1;
    scheme.score('A'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'A'_dna5) = -1;

    scheme.score('T'_dna5, 'C'_dna5) = -1;
    scheme.score('C'_dna5, 'T'_dna5) = -1;

    return scheme;
}

std::optional<CoOptimalPairwiseAligner::AlignmentResult> EvaluatedSplitRecords::getOptimalAlignment(
    const std::vector<CoOptimalPairwiseAligner::AlignmentResult> &alignResults,
    const double minComplementarity, const double minComplementarityFraction) {
    if (alignResults.empty()) {
        return std::nullopt;
    }

    std::optional<CoOptimalPairwiseAligner::AlignmentResult> optimalAlignment = std::nullopt;

    for (const auto &alignResult : alignResults) {
        if (alignResult.complementarity < minComplementarity ||
            alignResult.fraction < minComplementarityFraction) {
            continue;
        }

        if (!optimalAlignment.has_value() ||
            alignResult.complementarity > optimalAlignment.value().complementarity ||
            (alignResult.complementarity == optimalAlignment.value().complementarity &&
             alignResult.fraction > optimalAlignment.value().fraction)) {
            optimalAlignment = alignResult;
        }
    }

    return optimalAlignment;
}

std::optional<HybridizationResult> EvaluatedSplitRecords::calculateHybridization(
    std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
    const double mfeThreshold) {
    auto toString = [](const auto &seq) {
        return (seq | seqan3::views::to_char | seqan3::ranges::to<std::string>());
    };

    std::string interactionSeq = toString(seq1) + "&" + toString(seq2);

    vrna_fold_compound_t *vc =
        vrna_fold_compound(interactionSeq.c_str(), NULL, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
    char *structure = new char[interactionSeq.size() + 1];
    float mfe = vrna_cofold(interactionSeq.c_str(), structure);

    auto secondaryStructure = std::string(vrna_cut_point_insert(structure, seq1.size() + 1)) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    std::optional<CrosslinkingResult> crosslinkingResult =
        findCrosslinkingSites(seq1, seq2, secondaryStructure);

    delete[] structure;
    vrna_fold_compound_free(vc);

    if (mfe > mfeThreshold) {
        return std::nullopt;
    }

    return HybridizationResult{mfe, crosslinkingResult};
}

std::optional<CrosslinkingResult> EvaluatedSplitRecords::findCrosslinkingSites(
    std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
    std::vector<seqan3::dot_bracket3> &dotbracket) {
    if (seq1.empty() || seq2.empty() || dotbracket.empty()) {
        Logger::log(LogLevel::WARNING, "Empty input sequences or dot-bracket vector!");
        return std::nullopt;
    }

    std::vector<size_t> openPos;

    std::vector<NucleotidePairPositions> interactionPositions;
    interactionPositions.reserve(dotbracket.size() / 2);

    for (size_t i = 0; i < dotbracket.size(); ++i) {
        if (dotbracket[i].is_pair_open()) {
            openPos.push_back(i);
        } else if (dotbracket[i].is_pair_close()) {
            int open = openPos.back();
            openPos.pop_back();
            interactionPositions.emplace_back(open, i);
        }
    }

    // Sort pairing sites by first position
    std::sort(interactionPositions.begin(), interactionPositions.end(),
              [](const NucleotidePairPositions &a, const NucleotidePairPositions &b) {
                  return a.first < b.first;
              });

    if (!openPos.empty()) {
        Logger::log(LogLevel::WARNING, "Unpaired bases in dot-bracket notation! (" +
                                           std::to_string(openPos.size()) + ")");
        return std::nullopt;
    }

    if (interactionPositions.size() == 0) {
        return CrosslinkingResult{0, 0, 0, 0};
    }

    std::vector<NucleotidePositionsWindow> crosslinkingSites;
    crosslinkingSites.reserve(interactionPositions.size() - 1);

    size_t intraCrosslinkingScore = 0;
    size_t interCrosslinkingScore = 0;

    int preferredCrosslinkingCount = 0;
    int nonPreferredCrosslinkingCount = 0;
    int wobbleCrosslinkingCount = 0;

    // Iterate over all pairs of interaction sites until the second last pair (last pair
    // has no following pair for window)
    for (size_t i = 0; i < interactionPositions.size() - 1; ++i) {
        const NucleotidePositionsWindow window =
            std::make_pair(interactionPositions[i], interactionPositions[i + 1]);

        std::optional<InteractionWindow> interactionWindow =
            getContinuosNucleotideWindows(seq1, seq2, window);

        if (interactionWindow) {
            InteractionWindow &v_interactionWindow = interactionWindow.value();

            const NucleotideWindowPair nucleotidePairs =
                std::make_pair(v_interactionWindow.forwardWindowNucleotides,
                               v_interactionWindow.reverseWindowNucleotides);

            const auto it = crosslinkingScoringScheme.find(nucleotidePairs);
            if (it != crosslinkingScoringScheme.end()) {
                const size_t score = it->second;
                if (v_interactionWindow.isInterFragment) {
                    interCrosslinkingScore += score;
                    if (score == 3) {
                        preferredCrosslinkingCount++;
                    } else if (score == 2) {
                        nonPreferredCrosslinkingCount++;
                    } else if (score == 1) {
                        wobbleCrosslinkingCount++;
                    }
                } else {
                    intraCrosslinkingScore += score;
                }
                crosslinkingSites.emplace_back(window);
            }
        }
    }

    size_t minSeqLength = std::min(seq1.size(), seq2.size());

    if (minSeqLength == 0) {
        Logger::log(LogLevel::WARNING, "Division by zero!");
        return std::nullopt;
    }

    double normCrosslinkingScore = static_cast<double>(interCrosslinkingScore) / minSeqLength;
    return CrosslinkingResult{normCrosslinkingScore, preferredCrosslinkingCount,
                              nonPreferredCrosslinkingCount, wobbleCrosslinkingCount};
}

std::optional<InteractionWindow> EvaluatedSplitRecords::getContinuosNucleotideWindows(
    std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
    NucleotidePositionsWindow positionsPair) {
    std::pair<uint16_t, uint16_t> forwardPair =
        std::make_pair(positionsPair.first.first, positionsPair.second.first);
    std::pair<uint16_t, uint16_t> reversePair =
        std::make_pair(positionsPair.first.second, positionsPair.second.second);

    // Check that both pairs are from a continuous region
    if (forwardPair.first + 1 != forwardPair.second ||
        reversePair.first - 1 != reversePair.second) {
        return std::nullopt;
    }

    size_t seq1Length = seq1.size();

    bool forwardPairSplit = forwardPair.first < seq1Length && forwardPair.second >= seq1Length;
    bool reversePairSplit = reversePair.first >= seq1Length && reversePair.second < seq1Length;

    if (forwardPairSplit || reversePairSplit) {
        return std::nullopt;
    }

    bool firstPairInSeq1 = forwardPair.second < seq1Length;
    bool secondPairInSeq1 = reversePair.first < seq1Length;

    // Both pairs are in the first sequence
    if (firstPairInSeq1 && secondPairInSeq1) {
        return InteractionWindow{
            seqan3::dna5_vector{seq1[forwardPair.first], seq1[forwardPair.second]},
            seqan3::dna5_vector{seq1[reversePair.first], seq1[reversePair.second]}, forwardPair,
            reversePair, false};
    }

    // Both pairs are in the second sequence
    if (!firstPairInSeq1 && !secondPairInSeq1) {
        return InteractionWindow{seqan3::dna5_vector{seq2[forwardPair.first - seq1Length - 1],
                                                     seq2[forwardPair.second - seq1Length - 1]},
                                 seqan3::dna5_vector{seq2[reversePair.first - seq1Length - 1],
                                                     seq2[reversePair.second - seq1Length - 1]},
                                 forwardPair, reversePair, false};
    }

    // The first pair is in the first sequence and the second pair is in the second sequence
    if (firstPairInSeq1 && !secondPairInSeq1) {
        return InteractionWindow{
            seqan3::dna5_vector{seq1[forwardPair.first], seq1[forwardPair.second]},
            seqan3::dna5_vector{seq2[reversePair.first - seq1Length - 1],
                                seq2[reversePair.second - seq1Length - 1]},
            forwardPair, reversePair, true};
    }

    // Due to sorting of base pairs first pair in second sequence and second pair in first
    // sequence is not possible
    return std::nullopt;
}