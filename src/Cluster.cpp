#include "Cluster.hpp"

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
    return dtp::GenomicRegion{referenceIDs[referenceIDIndex], start, end};
}

dtp::Feature Segment::toFeature(const std::deque<std::string> &referenceIDs,
                                const std::string &featureID,
                                const std::string &featureType) const {
    return dtp::Feature{referenceIDs[referenceIDIndex], featureType, start, end, strand, featureID};
}

void Segment::merge(const Segment &other) {
    start = std::min(start, other.start);
    end = std::max(end, other.end);
    complementarityScore = std::max(complementarityScore, other.complementarityScore);
    hybridizationEnergy = std::min(hybridizationEnergy, other.hybridizationEnergy);
}

// ReadCluster
ReadCluster::ReadCluster(std::pair<Segment, Segment> segments,
                         std::vector<double> complementarityScores,
                         std::vector<double> hybridizationEnergies)
    : segments(segments),
      complementarityScores(complementarityScores),
      hybridizationEnergies(hybridizationEnergies) {}

std::optional<std::pair<Segment, Segment>> ReadCluster::getSortedElements(const Segment &segment1,
                                                                          const Segment &segment2) {
    if (segment1.recordID != segment2.recordID) {
        Logger::log(LogLevel::WARNING, "Record IDs do not match: ", segment1.recordID, " vs. ",
                    segment2.recordID, ". Make sure reads are sorted by name.");
        return std::nullopt;
    }
    return segment1.start < segment2.start ? std::make_pair(segment1, segment2)
                                           : std::make_pair(segment2, segment1);
}

std::optional<ReadCluster> ReadCluster::fromSegments(const Segment &segment1,
                                                     const Segment &segment2) {
    const auto sortedElements = getSortedElements(segment1, segment2);
    if (!sortedElements.has_value()) {
        return std::nullopt;
    }
    assert(sortedElements.value().first.complementarityScore ==
           sortedElements.value().second.complementarityScore);
    assert(sortedElements.value().first.hybridizationEnergy ==
           sortedElements.value().second.hybridizationEnergy);
    return ReadCluster{sortedElements.value(),
                       {sortedElements->first.complementarityScore},
                       {sortedElements->first.hybridizationEnergy}};
}

bool ReadCluster::operator<(const ReadCluster &a) const {
    return std::tie(segments.first.start, segments.second.start) <
           std::tie(a.segments.first.start, a.segments.second.start);
}

bool ReadCluster::overlaps(const ReadCluster &other, const int graceDistance) const {
    const bool firstOverlaps = segments.first.end + graceDistance >= other.segments.first.start &&
                               segments.first.start <= other.segments.first.end + graceDistance;

    const bool secondOverlaps =
        segments.second.end + graceDistance >= other.segments.second.start &&
        segments.second.start <= other.segments.second.end + graceDistance;

    return segments.first.referenceIDIndex == other.segments.first.referenceIDIndex &&
           segments.second.referenceIDIndex == other.segments.second.referenceIDIndex &&
           firstOverlaps && secondOverlaps &&
           segments.first.strand == other.segments.first.strand &&
           segments.second.strand == other.segments.second.strand;
}

void ReadCluster::merge(const ReadCluster &other) {
    segments.first.merge(other.segments.first);
    segments.second.merge(other.segments.second);
    complementarityScores.insert(complementarityScores.end(), other.complementarityScores.begin(),
                                 other.complementarityScores.end());
    hybridizationEnergies.insert(hybridizationEnergies.end(), other.hybridizationEnergies.begin(),
                                 other.hybridizationEnergies.end());
    count += other.count;
}

double ReadCluster::complementarityStatistics() const {
    assert(!complementarityScores.empty());
    assert(segments.first.complementarityScore == segments.second.complementarityScore);
    return helper::calculateMedian(complementarityScores) / segments.first.complementarityScore;
}

double ReadCluster::hybridizationEnergyStatistics() const {
    assert(!hybridizationEnergies.empty());
    assert(segments.first.hybridizationEnergy == segments.second.hybridizationEnergy);
    return std::sqrt(helper::calculateMedian(hybridizationEnergies) /
                     segments.first.hybridizationEnergy);
}

// Cluster
Cluster::Cluster(po::variables_map params)
    : params(params),
      featureAnnotator(params["features"].as<std::string>(),
                       params["featuretypes"].as<std::vector<std::string>>()) {}

void Cluster::iterateSplitsFile(const fs::path &splitsInPath,
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

    int subsetChunkSize = 10;
    std::vector<ReadCluster> subset;
    subset.reserve(subsetChunkSize);

    std::vector<ReadCluster> mergedClusters;

    for (auto &&records : splitsIn | seqan3::views::chunk(2)) {
        std::optional<Segment> segment1 = Segment::fromSamRecord(*records.begin());
        std::optional<Segment> segment2 = Segment::fromSamRecord(*(++records.begin()));

        if (!segment1 || !segment2) {
            continue;
        }

        auto cluster = ReadCluster::fromSegments(*segment1, *segment2);

        if (!cluster) {
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

    writeClustersToFile(mergedClusters, clusterOutPath, referenceIDs);
}

void Cluster::assignClustersToTranscripts(std::vector<ReadCluster> &clusters,
                                          const std::deque<std::string> &referenceIDs,
                                          const fs::path &unassignedSingletonsInPath,
                                          const fs::path &fragmentCountsInPath,
                                          const fs::path &supplementaryFeaturesOutPath,
                                          const fs::path &transcriptCountsOutPath) {
    std::ofstream supplementaryFeaturesOut(supplementaryFeaturesOutPath.string());

    if (!supplementaryFeaturesOut.is_open()) {
        Logger::log(LogLevel::ERROR,
                    "Could not open file: ", supplementaryFeaturesOutPath.string());
        exit(EXIT_FAILURE);
    }

    supplementaryFeaturesOut << "##gff-version 3\n";

    size_t supplementaryFeatureCount = 0;
    auto writeSupplementaryFeature = [&](const Segment &segment,
                                         std::stringstream &stream) -> std::string {
        const std::string transcriptID =
            "supplementary_" + std::to_string(supplementaryFeatureCount++);
        stream << referenceIDs[segment.referenceIDIndex] << "\t"
               << "RNAnue"
               << "\t"
               << "partial_transcript"
               << "\t" << segment.start << "\t" << segment.end << "\t"
               << "."
               << "\t" << static_cast<char>(segment.strand) << "\t"
               << "."
               << "\t"
               << "ID=" << transcriptID << ";"
               << "\n";
        return transcriptID;
    };

    std::unordered_map<std::string, size_t> transcriptCounts;
    dtp::FeatureMap supplementaryFeatureMap;

    std::stringstream supplementaryFeaturesStream;
    supplementaryFeaturesStream << "##gff-version 3\n";

    for (auto &cluster : clusters) {
        const auto &firstSegment = cluster.segments.first;
        const auto &secondSegment = cluster.segments.second;

        auto processSegment = [&](const Segment &segment) -> std::string {
            const auto segmentFeature =
                featureAnnotator.getBestOverlappingFeature(segment.toGenomicRegion(referenceIDs));
            if (!segmentFeature) {
                const std::string transcriptID =
                    writeSupplementaryFeature(segment, supplementaryFeaturesStream);
                transcriptCounts[transcriptID] = cluster.count;
                return transcriptID;
            } else {
                transcriptCounts[segmentFeature.value().id] += cluster.count;
                return segmentFeature.value().id;
            }
        };

        const std::string firstSegmentTranscriptID = processSegment(firstSegment);
        const std::string secondSegmentTranscriptID = processSegment(secondSegment);

        cluster.transcriptIDs = std::make_pair(firstSegmentTranscriptID, secondSegmentTranscriptID);
    }

    // Write supplementary features to file
    supplementaryFeaturesOut << supplementaryFeaturesStream.str();

    Annotation::FeatureAnnotator supplementaryFeatureAnnotator(supplementaryFeatureMap);

    assignNonAnnotatedSingletonsToSupplementaryFeatures(
        unassignedSingletonsInPath, supplementaryFeatureAnnotator, transcriptCounts);

    const size_t totalFragmentCount = parseSampleFragmentCount(fragmentCountsInPath);
    assignPValuesToClusters(clusters, transcriptCounts, totalFragmentCount);

    std::ofstream transcriptCountsOut(transcriptCountsOutPath.string());

    if (!transcriptCountsOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", transcriptCountsOutPath.string());
        exit(EXIT_FAILURE);
    }

    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptCountsOut << transcriptID << "\t" << count << "\n";
    }
}

std::unordered_map<std::string, double> Cluster::getTranscriptProbabilities(
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

void Cluster::assignPAdjustedValuesToClusters(std::vector<ReadCluster> &clusters) {
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

void Cluster::assignPValuesToClusters(
    std::vector<ReadCluster> &clusters,
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

void Cluster::assignNonAnnotatedSingletonsToSupplementaryFeatures(
    const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    seqan3::sam_file_input unassignedSingletonsIn{unassignedSingletonsInPath.string(),
                                                  sam_field_ids{}};

    for (auto &&records : unassignedSingletonsIn) {
        const auto region =
            dtp::GenomicRegion::fromSamRecord(records, unassignedSingletonsIn.header().ref_ids());

        if (!region) {
            continue;
        }

        const auto bestFeature = featureAnnotator.getBestOverlappingFeature(region.value());

        if (!bestFeature) {
            continue;
        }

        const std::string &transcriptID = bestFeature->id;
        assert(transcriptCounts.contains(transcriptID));
        ++transcriptCounts[transcriptID];
    }
}

void Cluster::writeClustersToFile(const std::vector<ReadCluster> &mergedClusters,
                                  const fs::path &clusterOutPath,
                                  const std::deque<std::string> &referenceIDs) {
    Logger::log(LogLevel::INFO, "Writing clusters to file: ", clusterOutPath.string());

    const auto getReferenceID = [&referenceIDs](const uint32_t referenceIDIndex) -> std::string {
        return referenceIDs[referenceIDIndex];
    };

    std::ofstream outputFile(clusterOutPath.string());
    if (!outputFile.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", clusterOutPath.string());
        exit(EXIT_FAILURE);
    }

    outputFile << "cluster_ID\tfst_seg_chr\tfst_seg_strd\tfst_seg_strt\tfst_seg_end\t"
                  "sec_seg_chr\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\tno_splits\t"
                  "fst_seg_len\tsec_seg_len\tgcs\tghs\tp_value\tpadj_value\n";

    size_t clusterID = 1;
    for (const auto &cluster : mergedClusters) {
        if (cluster.count < 4) {
            continue;
        }

        outputFile << "cluster" << clusterID++ << "\t";

        outputFile << getReferenceID(cluster.segments.first.referenceIDIndex) << "\t";
        outputFile << static_cast<char>(cluster.segments.first.strand) << "\t";
        outputFile << cluster.segments.first.start << "\t";
        outputFile << cluster.segments.first.end << "\t";

        outputFile << getReferenceID(cluster.segments.second.referenceIDIndex) << "\t";
        outputFile << static_cast<char>(cluster.segments.second.strand) << "\t";
        outputFile << cluster.segments.second.start << "\t";
        outputFile << cluster.segments.second.end << "\t";

        outputFile << cluster.count << "\t";
        outputFile << (cluster.segments.first.end + 1) - cluster.segments.first.start << "\t";
        outputFile << (cluster.segments.second.end + 1) - cluster.segments.second.start << "\t";
        outputFile << cluster.complementarityStatistics() << "\t";
        outputFile << cluster.hybridizationEnergyStatistics() << "\t";
        outputFile << cluster.pValue.value_or(1.0) << "\t";
        outputFile << cluster.pAdj.value_or(1.0);
        outputFile << "\n";
    }
}

void Cluster::mergeOverlappingClusters(std::vector<ReadCluster> &clusters) {
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

size_t Cluster::parseSampleFragmentCount(const fs::path &sampleCountsInPath) {
    std::ifstream sampleCountsIn(sampleCountsInPath);

    if (!sampleCountsIn.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", sampleCountsInPath.string());
        exit(EXIT_FAILURE);
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

void Cluster::start(pt::ptree sample) {
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