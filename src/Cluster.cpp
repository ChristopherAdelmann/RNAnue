#include "Cluster.hpp"

Cluster::Cluster(po::variables_map params)
    : params(params),
      featureAnnotator(params["features"].as<std::string>(),
                       params["featuretypes"].as<std::vector<std::string>>()) {}

void Cluster::iterate(const std::string &splitsInPath, const fs::path &unassignedSingletonsInPath,
                      const fs::path &clusterOutPath, const fs::path &supplementaryFeaturesOutPath,
                      const fs::path &clusterTranscriptCountsOutPath) {
    // input .sam file of sngl splits
    seqan3::sam_file_input splitsIn{splitsInPath, sam_field_ids{}};

    std::vector<std::optional<int32_t>> refOffsets;

    auto &header = splitsIn.header();
    std::vector<size_t> ref_lengths{};
    std::ranges::transform(header.ref_id_info, std::back_inserter(ref_lengths),
                           [](auto const &info) { return std::get<0>(info); });

    const std::deque<std::string> &referenceIDs = header.ref_ids();

    int subsetChunkSize = 5;
    std::vector<ReadCluster> subset;
    subset.reserve(subsetChunkSize);

    for (auto &&records : splitsIn | seqan3::views::chunk(2)) {
        std::optional<Segment> segment1 = Segment::fromSamRecord(*records.begin());
        std::optional<Segment> segment2 = Segment::fromSamRecord(*(++records.begin()));

        if (!segment1.has_value() || !segment2.has_value()) {
            continue;
        }

        auto cluster = ReadCluster::fromSegments(segment1.value(), segment2.value());

        if (!cluster.has_value()) {
            continue;
        }

        subset.push_back(std::move(cluster.value()));

        if (subset.size() == subsetChunkSize) {
            mergeOverlappingClusters(subset);
            result.insert(result.end(), std::make_move_iterator(subset.begin()),
                          std::make_move_iterator(subset.end()));
            subset.clear();
        }
    }

    if (!subset.empty()) {
        mergeOverlappingClusters(subset);
        result.insert(result.end(), std::make_move_iterator(subset.begin()),
                      std::make_move_iterator(subset.end()));
    }

    mergeOverlappingClusters(result);

    assignClustersToTranscripts(result, referenceIDs, unassignedSingletonsInPath,
                                supplementaryFeaturesOutPath, clusterTranscriptCountsOutPath);

    writeClustersToFile(clusterOutPath, referenceIDs);
}

void Cluster::assignClustersToTranscripts(const std::vector<ReadCluster> &clusters,
                                          const std::deque<std::string> &referenceIDs,
                                          const fs::path &unassignedSingletonsInPath,
                                          const fs::path &supplementaryFeaturesOutPath,
                                          const fs::path &transcriptCountsOutPath) {
    std::ofstream supplementaryFeaturesOut(supplementaryFeaturesOutPath.string());

    if (!supplementaryFeaturesOut.is_open()) {
        Logger::log(LogLevel::ERROR,
                    "Could not open file: ", supplementaryFeaturesOutPath.string());
        exit(1);
    }

    supplementaryFeaturesOut << "##gff-version 3\n";

    size_t supplementaryFeatureCount = 0;
    auto writeSupplementaryFeature = [&](const Segment &segment) -> std::string {
        const std::string transcriptID =
            "supplementary_" + std::to_string(supplementaryFeatureCount++);
        supplementaryFeaturesOut << referenceIDs[segment.referenceIDIndex] << "\t";
        supplementaryFeaturesOut << "Cluster"
                                 << "\t";
        supplementaryFeaturesOut << "partial_transcript"
                                 << "\t";
        supplementaryFeaturesOut << segment.start << "\t";
        supplementaryFeaturesOut << segment.end << "\t";
        supplementaryFeaturesOut << "."
                                 << "\t";
        supplementaryFeaturesOut << static_cast<char>(segment.strand) << "\t";
        supplementaryFeaturesOut << "."
                                 << "\t";
        supplementaryFeaturesOut << "ID=" << transcriptID << ";"
                                 << "\n";
        return transcriptID;
    };

    std::unordered_map<std::string, size_t> transcriptCounts;
    std::vector<std::pair<std::string, std::string>> unassignedSingletonsInPath;

    dtp::FeatureMap supplementaryFeatureMap;

    auto addSegmentToSupplementaryFeatures = [&](const Segment &segment, const size_t count) {
        const std::string transcriptID = writeSupplementaryFeature(segment);
        const dtp::Feature feature = segment.toFeature(referenceIDs, transcriptID);
        supplementaryFeatureMap[feature.referenceID].push_back(feature);
        transcriptCounts[transcriptID] = count;
        return transcriptID;
    };

    for (auto &cluster : clusters) {
        const auto &firstSegment = cluster.segments.first;
        const auto &secondSegment = cluster.segments.second;

        const auto firstSegmentFeature =
            featureAnnotator.getBestOverlappingFeature(firstSegment.toGenomicRegion(referenceIDs));
        const auto secondSegmentFeature =
            featureAnnotator.getBestOverlappingFeature(secondSegment.toGenomicRegion(referenceIDs));

        // TODO Segments from different clusters that overlap should be merged to one supplementary
        // feature
        std::string firstSegmentTranscriptID;
        if (!firstSegmentFeature.has_value()) {
            firstSegmentTranscriptID =
                addSegmentToSupplementaryFeatures(firstSegment, cluster.count);
        } else {
            transcriptCounts[firstSegmentFeature.value().id] += cluster.count;
            firstSegmentTranscriptID = firstSegmentFeature.value().id;
        }

        std::string secondSegmentTranscriptID;
        if (!secondSegmentFeature.has_value()) {
            secondSegmentTranscriptID =
                addSegmentToSupplementaryFeatures(secondSegment, cluster.count);
        } else {
            transcriptCounts[secondSegmentFeature.value().id] += cluster.count;
            secondSegmentTranscriptID = secondSegmentFeature.value().id;
        }
    }

    Annotation::FeatureAnnotator supplementaryFeatureAnnotator(supplementaryFeatureMap);

    assignUnassignedSingletonsToSupplementaryFeatures(
        unassignedSingletonsInPath, supplementaryFeatureAnnotator, transcriptCounts);

    std::ofstream transcriptCountsOut(transcriptCountsOutPath.string());

    if (!transcriptCountsOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", transcriptCountsOutPath.string());
        exit(1);
    }

    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptCountsOut << transcriptID << "\t" << count << "\n";
    }
}

void Cluster::assignUnassignedSingletonsToSupplementaryFeatures(
    const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    seqan3::sam_file_input unassignedSingletonsIn{unassignedSingletonsInPath.string(),
                                                  sam_field_ids{}};

    for (auto &&records : unassignedSingletonsIn) {
        dtp::GenomicRegion region =
            dtp::GenomicRegion::fromSamRecord(records, unassignedSingletonsIn.header().ref_ids());

        const auto bestFeature = featureAnnotator.getBestOverlappingFeature(region);

        if (!bestFeature.has_value()) {
            continue;
        }

        const std::string &transcriptID = bestFeature.value().id;
        assert(transcriptCounts.contains(transcriptID));
        ++transcriptCounts[transcriptID];
    }
}

void Cluster::writeClustersToFile(const fs::path &clusterOutPath,
                                  const std::deque<std::string> &referenceIDs) {
    Logger::log(LogLevel::INFO, "Writing clusters to file: ", clusterOutPath.string());

    const auto getReferenceID = [&referenceIDs](const uint32_t referenceIDIndex) -> std::string {
        return referenceIDs[referenceIDIndex];
    };

    std::ofstream outputFile(clusterOutPath.string());
    if (!outputFile.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", clusterOutPath.string());
        exit(1);
    }

    outputFile << "cluster_ID\tfst_seg_chr\tfst_seg_strd\tfst_seg_strt\tfst_seg_end\t"
                  "sec_seg_chr\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\tno_splits\t"
                  "fst_seg_len\tsec_seg_len\tgcs\tghs\tcmplScores\n";

    size_t clusterID = 1;
    for (const auto &cluster : result) {
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
        for (const auto complementarityScore : cluster.complementarityScores) {
            outputFile << complementarityScore << ",";
        }
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

void Cluster::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    pt::ptree output = sample.get_child("output");

    std::string splitsInPath = input.get<std::string>("splits");
    fs::path unassignedSingletonsInPath = input.get<std::string>("singletonunassigned");
    fs::path clusterOutPath = output.get<std::string>("clusters");
    fs::path supplementaryFeaturesOutPath = output.get<std::string>("supplementaryfeatures");
    fs::path clusterTranscriptCountsOutPath = output.get<std::string>("clustertranscriptcounts");

    iterate(splitsInPath, unassignedSingletonsInPath, clusterOutPath, supplementaryFeaturesOutPath,
            clusterTranscriptCountsOutPath);
}
