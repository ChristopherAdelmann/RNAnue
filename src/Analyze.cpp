#include "Analyze.hpp"

// Segment
std::optional<Segment> Segment::fromSamRecord(const dtp::SamRecord &record) {
    if (!record.reference_position().has_value() || !record.reference_id().has_value()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const dtp::Strand strand{isReverseStrand ? dtp::Strand::REVERSE : dtp::Strand::FORWARD};

    const auto start = record.reference_position();
    const std::optional<int32_t> end = recordEndPosition(record);

    if (!start.has_value() || !end.has_value()) {
        return std::nullopt;
    }

    const double hybridizationEnergy = record.tags().get<"XE"_tag>();
    const double complementarityScore = record.tags().get<"XC"_tag>();

    return Segment{record.id(),
                   record.reference_id().value(),
                   strand,
                   start.value(),
                   end.value(),
                   complementarityScore,
                   hybridizationEnergy};
}

dtp::GenomicRegion Segment::toGenomicRegion(const std::deque<std::string> &referenceIDs) const {
    return dtp::GenomicRegion{referenceIDs[referenceIDIndex], start, end, strand};
}

dtp::Feature Segment::toFeature(const std::deque<std::string> &referenceIDs,
                                const std::string &featureID,
                                const std::string &featureType) const {
    return dtp::Feature{referenceIDs[referenceIDIndex], featureType, start, end, strand, featureID};
}

void Segment::merge(const Segment &other) {
    start = std::min(start, other.start);
    end = std::max(end, other.end);
    maxComplementarityScore = std::max(maxComplementarityScore, other.maxComplementarityScore);
    minHybridizationEnergy = std::min(minHybridizationEnergy, other.minHybridizationEnergy);
}

// InteractionCluster
InteractionCluster::InteractionCluster(std::pair<Segment, Segment> segments,
                                       const std::vector<double> &complementarityScores,
                                       const std::vector<double> &hybridizationEnergies)
    : segments(segments),
      complementarityScores(complementarityScores),
      hybridizationEnergies(hybridizationEnergies) {}

std::optional<std::pair<Segment, Segment>> InteractionCluster::getSortedElements(
    const Segment &segment1, const Segment &segment2) {
    if (segment1.recordID != segment2.recordID) {
        Logger::log(LogLevel::WARNING, "Record IDs do not match: ", segment1.recordID, " vs. ",
                    segment2.recordID, ". Make sure reads are sorted by name.");
        return std::nullopt;
    }

    // Optimized version
    return segment1.referenceIDIndex < segment2.referenceIDIndex ||
                   (segment1.referenceIDIndex == segment2.referenceIDIndex &&
                    segment1.start < segment2.start)
               ? std::make_pair(segment1, segment2)
               : std::make_pair(segment2, segment1);
}

std::optional<InteractionCluster> InteractionCluster::fromSegments(const Segment &segment1,
                                                                   const Segment &segment2) {
    const auto sortedElements = getSortedElements(segment1, segment2);
    if (!sortedElements.has_value()) [[unlikely]] {
        return std::nullopt;
    }

    assert(sortedElements.value().first.maxComplementarityScore ==
           sortedElements.value().second.maxComplementarityScore);
    assert(sortedElements.value().first.minHybridizationEnergy ==
           sortedElements.value().second.minHybridizationEnergy);

    return InteractionCluster{sortedElements.value(),
                              {sortedElements->first.maxComplementarityScore},
                              {sortedElements->first.minHybridizationEnergy}};
}

bool InteractionCluster::operator<(const InteractionCluster &a) const {
    return std::tie(segments.first.start, segments.second.start) <
           std::tie(a.segments.first.start, a.segments.second.start);
}

bool InteractionCluster::overlaps(const InteractionCluster &other, const int graceDistance) const {
    bool isSameReferenceAndStrand =
        segments.first.referenceIDIndex == other.segments.first.referenceIDIndex &&
        segments.second.referenceIDIndex == other.segments.second.referenceIDIndex &&
        segments.first.strand == other.segments.first.strand &&
        segments.second.strand == other.segments.second.strand;

    if (!isSameReferenceAndStrand) {
        return false;
    }

    const bool firstOverlaps = segments.first.end + graceDistance >= other.segments.first.start &&
                               segments.first.start <= other.segments.first.end + graceDistance;
    const bool secondOverlaps =
        segments.second.end + graceDistance >= other.segments.second.start &&
        segments.second.start <= other.segments.second.end + graceDistance;

    return firstOverlaps && secondOverlaps;
}

void InteractionCluster::merge(const InteractionCluster &other) {
    assert(segments.first.referenceIDIndex == other.segments.first.referenceIDIndex);
    assert(segments.second.referenceIDIndex == other.segments.second.referenceIDIndex);
    assert(segments.first.strand == other.segments.first.strand);
    assert(segments.second.strand == other.segments.second.strand);

    segments.first.merge(other.segments.first);
    segments.second.merge(other.segments.second);

    complementarityScores.reserve(complementarityScores.size() +
                                  other.complementarityScores.size());
    hybridizationEnergies.reserve(hybridizationEnergies.size() +
                                  other.hybridizationEnergies.size());

    complementarityScores.insert(complementarityScores.end(),
                                 std::make_move_iterator(other.complementarityScores.begin()),
                                 std::make_move_iterator(other.complementarityScores.end()));
    hybridizationEnergies.insert(hybridizationEnergies.end(),
                                 std::make_move_iterator(other.hybridizationEnergies.begin()),
                                 std::make_move_iterator(other.hybridizationEnergies.end()));

    count += other.count;
}

double InteractionCluster::complementarityStatistics() const {
    assert(!complementarityScores.empty());
    assert(segments.first.maxComplementarityScore == segments.second.maxComplementarityScore);
    return helper::calculateMedian(complementarityScores) * segments.first.maxComplementarityScore;
}

double InteractionCluster::hybridizationEnergyStatistics() const {
    assert(!hybridizationEnergies.empty());
    assert(segments.first.minHybridizationEnergy == segments.second.minHybridizationEnergy);
    return std::sqrt(helper::calculateMedian(hybridizationEnergies) *
                     segments.first.minHybridizationEnergy);
}

// Analyze
Analyze::Analyze(po::variables_map params)
    : params(params),
      featureAnnotator(params["features"].as<std::string>(),
                       params["featuretypes"].as<std::vector<std::string>>()) {}

void Analyze::iterateSplitsFile(const fs::path &splitsInPath,
                                const fs::path &unassignedSingletonsInPath,
                                const fs::path &fragmentCountsInPath,
                                const fs::path &clusterOutPath,
                                const fs::path &supplementaryFeaturesOutPath,
                                const fs::path &clusterTranscriptCountsOutPath) {
    seqan3::sam_file_input splitsIn{splitsInPath.string(), sam_field_ids{}};

    auto &header = splitsIn.header();
    std::vector<size_t> ref_lengths{};
    std::ranges::transform(header.ref_id_info, std::back_inserter(ref_lengths),
                           [](auto const &info) { return std::get<0>(info); });

    const std::deque<std::string> &referenceIDs = header.ref_ids();

    size_t subsetChunkSize = 10;
    std::vector<InteractionCluster> subset;
    subset.reserve(subsetChunkSize);

    std::vector<InteractionCluster> mergedClusters;

    for (auto &&records : splitsIn | seqan3::views::chunk(2)) {
        std::optional<Segment> segment1 = Segment::fromSamRecord(*records.begin());
        std::optional<Segment> segment2 = Segment::fromSamRecord(*(++records.begin()));

        if (!segment1 || !segment2) [[unlikely]] {
            continue;
        }

        auto cluster = InteractionCluster::fromSegments(*segment1, *segment2);

        if (!cluster) [[unlikely]] {
            continue;
        }

        subset.push_back(std::move(*cluster));

        if (subset.size() == subsetChunkSize) {
            mergeOverlappingClusters(subset);

            mergedClusters.reserve(mergedClusters.size() + subset.size());
            std::move(subset.begin(), subset.end(), std::back_inserter(mergedClusters));
            subset.clear();
        }
    }

    if (!subset.empty()) {
        mergeOverlappingClusters(subset);

        mergedClusters.reserve(mergedClusters.size() + subset.size());
        std::move(subset.begin(), subset.end(), std::back_inserter(mergedClusters));
    }

    mergeOverlappingClusters(mergedClusters);

    assignClustersToTranscripts(mergedClusters, referenceIDs, unassignedSingletonsInPath,
                                fragmentCountsInPath, supplementaryFeaturesOutPath,
                                clusterTranscriptCountsOutPath);

    writeInteractionsToFile(mergedClusters, clusterOutPath, referenceIDs);
}

void Analyze::assignClustersToTranscripts(std::vector<InteractionCluster> &clusters,
                                          const std::deque<std::string> &referenceIDs,
                                          const fs::path &unassignedSingletonsInPath,
                                          const fs::path &fragmentCountsInPath,
                                          const fs::path &supplementaryFeaturesOutPath,
                                          const fs::path &transcriptCountsOutPath) {
    std::unordered_map<std::string, size_t> transcriptCounts;
    std::unordered_map<std::string, std::string> mergedTranscriptIDsIntoTranscriptIDs;

    Annotation::FeatureAnnotator supplementaryFeatureAnnotator;

    const int mergeGraceDistance = params["clustdist"].as<int>();
    const auto annotationOrientation = params["orientation"].as<Annotation::Orientation>();

    for (auto &cluster : clusters) {
        const auto &firstSegment = cluster.segments.first;
        const auto &secondSegment = cluster.segments.second;

        auto processSegment = [&](const Segment &segment) -> std::string {
            const auto region = segment.toGenomicRegion(referenceIDs);

            const auto segmentFeature =
                featureAnnotator.getBestOverlappingFeature(region, annotationOrientation);

            if (segmentFeature.has_value()) {
                transcriptCounts[segmentFeature.value().id] += cluster.count;
                return segmentFeature.value().id;
            }

            const auto supplementarySegmentFeature =
                supplementaryFeatureAnnotator.getBestOverlappingFeature(
                    segment.toGenomicRegion(referenceIDs), annotationOrientation);

            if (supplementarySegmentFeature.has_value()) {
                transcriptCounts[supplementarySegmentFeature.value().id] += cluster.count;
                return supplementarySegmentFeature.value().id;
            }

            const auto mergeResult = supplementaryFeatureAnnotator.mergeInsert(
                segment.toGenomicRegion(referenceIDs), mergeGraceDistance);

            transcriptCounts[mergeResult.featureID] += cluster.count;

            for (const auto &mergedFeatureID : mergeResult.mergedFeatureIDs) {
                mergedTranscriptIDsIntoTranscriptIDs[mergedFeatureID] = mergeResult.featureID;
                transcriptCounts[mergeResult.featureID] += transcriptCounts[mergedFeatureID];
                transcriptCounts.erase(mergedFeatureID);
            }

            return mergeResult.featureID;
        };

        const std::string firstSegmentTranscriptID = processSegment(firstSegment);
        const std::string secondSegmentTranscriptID = processSegment(secondSegment);

        cluster.transcriptIDs = std::make_pair(firstSegmentTranscriptID, secondSegmentTranscriptID);
    }

    for (auto &cluster : clusters) {
        if (cluster.transcriptIDs.has_value()) {
            const auto &firstTranscriptID = cluster.transcriptIDs->first;
            const auto &secondTranscriptID = cluster.transcriptIDs->second;

            if (mergedTranscriptIDsIntoTranscriptIDs.contains(firstTranscriptID)) {
                cluster.transcriptIDs->first =
                    mergedTranscriptIDsIntoTranscriptIDs[firstTranscriptID];
            }

            if (mergedTranscriptIDsIntoTranscriptIDs.contains(secondTranscriptID)) {
                cluster.transcriptIDs->second =
                    mergedTranscriptIDsIntoTranscriptIDs[secondTranscriptID];
            }
        }
    }

    FeatureWriter::write(supplementaryFeatureAnnotator.getFeatureTreeMap(),
                         supplementaryFeaturesOutPath.string(), Annotation::FileType::GFF);

    assignNonAnnotatedSingletonsToSupplementaryFeatures(
        unassignedSingletonsInPath, supplementaryFeatureAnnotator, transcriptCounts);

    const size_t totalFragmentCount = parseSampleFragmentCount(fragmentCountsInPath);
    assignPValuesToClusters(clusters, transcriptCounts, totalFragmentCount);

    std::ofstream transcriptCountsOut(transcriptCountsOutPath.string());

    if (!transcriptCountsOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", transcriptCountsOutPath.string());
    }

    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptCountsOut << transcriptID << "\t" << count << "\n";
    }
}

std::unordered_map<std::string, double> Analyze::getTranscriptProbabilities(
    const std::unordered_map<std::string, size_t> &transcriptCounts,
    const size_t totalFragmentCount) {
    std::unordered_map<std::string, double> transcriptProbabilities;

    double cumulativeProbability = 0;
    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptProbabilities[transcriptID] = static_cast<double>(count) / totalFragmentCount;
        cumulativeProbability += transcriptProbabilities[transcriptID];
    }

    // Normalize probabilities to sum to 1
    for (auto &[transcriptID, probability] : transcriptProbabilities) {
        probability /= cumulativeProbability;
    }

    return transcriptProbabilities;
}

void Analyze::assignPAdjustedValuesToClusters(std::vector<InteractionCluster> &clusters) {
    std::vector<std::pair<double, size_t>> pValuesWithIndex;
    pValuesWithIndex.reserve(clusters.size());  // Reserve capacity to avoid reallocations

    for (size_t i = 0; i < clusters.size(); ++i) {
        if (clusters[i].pValue.has_value()) {
            pValuesWithIndex.emplace_back(clusters[i].pValue.value(), i);
        }
    }

    std::sort(pValuesWithIndex.begin(), pValuesWithIndex.end());

    const size_t numPValues = pValuesWithIndex.size();
    double minAdjPValue = 1.0;
    for (size_t i = numPValues; i-- > 0;) {
        double rank = static_cast<double>(i + 1);
        double adjustedPValue =
            std::min(minAdjPValue, (pValuesWithIndex[i].first * numPValues / rank));
        minAdjPValue = adjustedPValue;
        clusters[pValuesWithIndex[i].second].pAdj = adjustedPValue;
    }
}

void Analyze::assignPValuesToClusters(
    std::vector<InteractionCluster> &clusters,
    const std::unordered_map<std::string, size_t> &transcriptCounts,
    const size_t totalFragmentCount) {
    const auto transcriptProbabilities =
        getTranscriptProbabilities(transcriptCounts, totalFragmentCount);

    for (auto &cluster : clusters) {
        if (!cluster.transcriptIDs.has_value()) {
            Logger::log(LogLevel::WARNING, "Some cluster is missing transcript IDs.");
            continue;
        }

        const auto &firstTranscriptID = cluster.transcriptIDs->first;
        const auto &secondTranscriptID = cluster.transcriptIDs->second;

        auto findProbability = [&](const std::string &transcriptID) -> std::optional<double> {
            auto it = transcriptProbabilities.find(transcriptID);
            if (it != transcriptProbabilities.end()) {
                return it->second;
            }
            return std::nullopt;
        };

        auto firstTranscriptProbability = findProbability(firstTranscriptID);
        auto secondTranscriptProbability = findProbability(secondTranscriptID);

        if (!firstTranscriptProbability || !secondTranscriptProbability) {
            Logger::log(LogLevel::WARNING, "Could not find transcript probabilities for cluster.");
            continue;
        }

        const double ligationByChanceProbability =
            (firstTranscriptID == secondTranscriptID)
                ? firstTranscriptProbability.value() * secondTranscriptProbability.value()
                : 2 * firstTranscriptProbability.value() * secondTranscriptProbability.value();
        ;

        const auto binomialDistribution =
            math::binomial_distribution((double)totalFragmentCount, ligationByChanceProbability);

        const double pValue = 1 - math::cdf(binomialDistribution, cluster.count);

        cluster.pValue = pValue;
    }

    assignPAdjustedValuesToClusters(clusters);
}

void Analyze::assignNonAnnotatedSingletonsToSupplementaryFeatures(
    const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    seqan3::sam_file_input unassignedSingletonsIn{unassignedSingletonsInPath.string(),
                                                  sam_field_ids{}};

    const auto annotationOrientation = params["orientation"].as<Annotation::Orientation>();

    for (auto &&records : unassignedSingletonsIn) {
        const auto region =
            dtp::GenomicRegion::fromSamRecord(records, unassignedSingletonsIn.header().ref_ids());

        if (!region) {
            continue;
        }

        const auto bestFeature =
            featureAnnotator.getBestOverlappingFeature(region.value(), annotationOrientation);

        if (!bestFeature) {
            continue;
        }

        const std::string &transcriptID = bestFeature->id;
        assert(transcriptCounts.contains(transcriptID));
        ++transcriptCounts[transcriptID];
    }
}

void Analyze::writeBEDLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                                 const std::deque<std::string> &referenceIDs,
                                 std::ofstream &bedOut) {
    const std::string randomColor = helper::generateRandomHexColor();
    bedOut << getReferenceID(cluster.segments.first.referenceIDIndex, referenceIDs) << "\t";
    bedOut << cluster.segments.first.start << "\t";
    bedOut << cluster.segments.first.end + 1 << "\t";
    bedOut << "cluster" << clusterID << "_segment1"
           << "\t";
    bedOut << "0\t";
    bedOut << static_cast<char>(cluster.segments.first.strand) << "\t";
    bedOut << cluster.segments.first.start << "\t";
    bedOut << cluster.segments.first.end + 1 << "\t";
    bedOut << randomColor << "\n";

    bedOut << getReferenceID(cluster.segments.second.referenceIDIndex, referenceIDs) << "\t";
    bedOut << cluster.segments.second.start << "\t";
    bedOut << cluster.segments.second.end + 1 << "\t";
    bedOut << "cluster" << clusterID << "_segment2"
           << "\t";
    bedOut << "0\t";
    bedOut << static_cast<char>(cluster.segments.second.strand) << "\t";
    bedOut << cluster.segments.second.start << "\t";
    bedOut << cluster.segments.second.end + 1 << "\t";
    bedOut << randomColor << "\n";
}

void Analyze::writeInteractionLineToFile(const InteractionCluster &cluster,
                                         const std::string &clusterID,
                                         const std::deque<std::string> &referenceIDs,
                                         std::ofstream &interactionOut) {
    interactionOut << "cluster" << clusterID << "\t";

    interactionOut << cluster.transcriptIDs->first << "\t";
    interactionOut << getReferenceID(cluster.segments.first.referenceIDIndex, referenceIDs) << "\t";
    interactionOut << static_cast<char>(cluster.segments.first.strand) << "\t";
    interactionOut << cluster.segments.first.start << "\t";
    interactionOut << cluster.segments.first.end << "\t";

    interactionOut << cluster.transcriptIDs->second << "\t";
    interactionOut << getReferenceID(cluster.segments.second.referenceIDIndex, referenceIDs)
                   << "\t";
    interactionOut << static_cast<char>(cluster.segments.second.strand) << "\t";
    interactionOut << cluster.segments.second.start << "\t";
    interactionOut << cluster.segments.second.end << "\t";

    interactionOut << cluster.count << "\t";
    interactionOut << cluster.complementarityStatistics() << "\t";
    interactionOut << cluster.hybridizationEnergyStatistics() << "\t";
    interactionOut << cluster.pValue.value_or(1.0) << "\t";
    interactionOut << cluster.pAdj.value_or(1.0);
    interactionOut << "\n";
}

void Analyze::writeInteractionsToFile(const std::vector<InteractionCluster> &mergedClusters,
                                      const fs::path &clusterOutPath,
                                      const std::deque<std::string> &referenceIDs) {
    Logger::log(LogLevel::INFO, "Writing interactions to file: ", clusterOutPath.string());

    std::ofstream interactionOut(clusterOutPath.string());
    if (!interactionOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", clusterOutPath.string());
    }

    interactionOut << "cluster_ID\tfst_feat_id\tfst_seg_chr\tfst_seg_strd\tfst_seg_strt\tfst_seg_"
                      "end\tsec_feat_id\t"
                      "sec_seg_chr\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\tno_splits\t"
                      "gcs\tghs\tp_value\tpadj_value\n";

    fs::path bedOutPath = clusterOutPath;
    bedOutPath.replace_extension(".bed");

    std::ofstream bedOut(bedOutPath.string());
    if (!bedOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", bedOutPath.string());
    }

    const std::string sampleName = clusterOutPath.stem().string();
    bedOut << "track name=\"" << sampleName
           << " RNA-RNA interactions\" description=\"Segments of interacting RNA "
              "clusters derived from DDD-Experiment\" itemRgb=\"On\"\n";

    const double pValueThreshold = params["padj"].as<double>();
    const size_t clusterCountThreshold = params["mincount"].as<int>();

    auto isClusterValid = [&](const InteractionCluster &cluster) {
        return cluster.pAdj.has_value() && cluster.pAdj.value() <= pValueThreshold &&
               cluster.count >= clusterCountThreshold;
    };

    size_t clusterID = 0;
    for (const auto &cluster : mergedClusters) {
        if (!isClusterValid(cluster)) {
            continue;
        }
        writeBEDLineToFile(cluster, std::to_string(clusterID), referenceIDs, bedOut);
        writeInteractionLineToFile(cluster, std::to_string(clusterID), referenceIDs,
                                   interactionOut);
        ++clusterID;
    }
}

void Analyze::mergeOverlappingClusters(std::vector<InteractionCluster> &clusters) {
    // sort by first segment (and the second)
    std::ranges::sort(clusters, std::less{});

    int clusterDistance = params["clustdist"].as<int>();

    for (size_t i = 0; i < clusters.size(); ++i) {
        for (size_t j = i + 1; j < clusters.size(); ++j) {
            // checks if the first overlaps
            if (clusters[i].overlaps(clusters[j], clusterDistance)) {
                clusters[j].merge(clusters[i]);
                clusters.erase(clusters.begin() + i);  // remove cluster
                --i;
                break;
            }
        }
    }
}

size_t Analyze::parseSampleFragmentCount(const fs::path &sampleCountsInPath) {
    std::ifstream sampleCountsIn(sampleCountsInPath);

    if (!sampleCountsIn.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", sampleCountsInPath.string());
    }

    sampleCountsIn.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip header

    size_t totalTranscriptCount = 0;
    std::string line;
    while (std::getline(sampleCountsIn, line)) {
        std::istringstream iss(line);

        size_t column = 0;
        for (std::string token; std::getline(iss, token, '\t');) {
            if (column != 0) {
                totalTranscriptCount += std::stoul(token);
            }
            ++column;
        }

        break;  // Should only be one line
    }

    return totalTranscriptCount;
}

void Analyze::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    fs::path splitsInPath = input.get<std::string>("splits");
    fs::path unassignedSingletonsInPath = input.get<std::string>("singletonunassigned");
    fs::path fragmentCountsInPath = input.get<std::string>("samplecounts");

    fs::path clusterOutPath = output.get<std::string>("clusters");
    fs::path supplementaryFeaturesOutPath = output.get<std::string>("supplementaryfeatures");
    fs::path clusterTranscriptCountsOutPath = output.get<std::string>("clustertranscriptcounts");

    iterateSplitsFile(splitsInPath, unassignedSingletonsInPath, fragmentCountsInPath,
                      clusterOutPath, supplementaryFeaturesOutPath, clusterTranscriptCountsOutPath);
}
