//
// Created by Richard Albin Schaefer on 3/4/24.
//
#include "Cluster.hpp"

Cluster::Cluster(po::variables_map params) : params(params) {}

void Cluster::iterate(std::string splitsInPath) {
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

    writeClustersToFile(referenceIDs);
}

void Cluster::writeClustersToFile(const std::deque<std::string> &referenceIDs) {
    // retrieve output directory
    fs::path output = fs::path(params["outdir"].as<std::string>()) / constants::pipelines::CLUSTER;
    fs::path cluster_results = output / "clusters.tab";

    Logger::log(LogLevel::INFO, "Writing clusters to file: ", cluster_results.string());

    const auto getStrand = [](const uint32_t flag) -> std::string { return flag == 0 ? "+" : "-"; };
    const auto getReferenceID = [&referenceIDs](const uint32_t referenceIDIndex) -> std::string {
        return referenceIDs[referenceIDIndex];
    };

    std::ofstream outputFile(cluster_results.string());
    if (!outputFile.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", cluster_results.string());
        exit(1);
    }

    outputFile << "clustID\tfst_seg_chr\tfst_seg_strd\tfst_seg_strt\tfst_seg_end\t"
                  "sec_seg_chr\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\tno_splits\t"
                  "fst_seg_len\tsec_seg_len\n";

    size_t clusterID = 1;
    for (const auto &cluster : result) {
        outputFile << "cluster" << clusterID++ << "\t";

        outputFile << getReferenceID(cluster.elements.first.referenceIDIndex) << "\t";
        outputFile << getStrand(cluster.elements.first.flag) << "\t";
        outputFile << cluster.elements.first.start << "\t";
        outputFile << cluster.elements.first.end << "\t";

        outputFile << getReferenceID(cluster.elements.second.referenceIDIndex) << "\t";
        outputFile << getStrand(cluster.elements.second.flag) << "\t";
        outputFile << cluster.elements.second.start << "\t";
        outputFile << cluster.elements.second.end << "\t";

        outputFile << cluster.count << "\t";
        outputFile << (cluster.elements.first.end + 1) - cluster.elements.first.start << "\t";
        outputFile << (cluster.elements.second.end + 1) - cluster.elements.second.start << "\n";
    }

    outputFile.close();
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

    std::string splits = input.get<std::string>("splits");

    iterate(splits);
}
