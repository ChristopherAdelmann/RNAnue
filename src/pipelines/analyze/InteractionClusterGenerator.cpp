#include "InteractionClusterGenerator.hpp"

namespace pipelines {
namespace analyze {

/**
 * @brief Merges overlapping interaction clusters and returns these.
 *
 * This function takes a list of interaction clusters, sorts them from back to front,
 * and merges any overlapping clusters. This is done from back to front while closing clusters that
 * are further back than the current cluster.
 *
 * @param clusters A reference to the list of interaction clusters to be merged.
 * @return A list of finalized interaction clusters.
 */
InteractionClusterGenerator::InteractionClusters InteractionClusterGenerator::mergeClusters(
    InteractionClusterGenerator::InteractionClusters& clusters) noexcept {
    // Clusters should be sorted from back to front
    std::sort(std::execution::par_unseq, clusters.begin(), clusters.end(), std::less<>{});

    for (auto& cluster : clusters | std::views::reverse) {
        if (openClusterQueue.empty()) {
            openClusterQueue.emplace_front(std::move(cluster));
            continue;
        }

        bool clusterWasMerged = false;

        auto prevIter = openClusterQueue.before_begin();

        for (auto iter = openClusterQueue.begin(); iter != openClusterQueue.end();) {
            if (iter->overlaps(cluster, graceDistance)) {
                iter->merge(cluster);
                clusterWasMerged = true;
                break;
            } else if (clusterIsBeforeOpenCluster(cluster, *iter)) {
                finalizeCluster(std::move(*iter));
                iter = openClusterQueue.erase_after(prevIter);
                continue;
            }

            prevIter = iter;
            ++iter;
        }

        if (!clusterWasMerged) {
            openClusterQueue.emplace_after(prevIter, std::move(cluster));
        }
    }

    for (auto& cluster : openClusterQueue) {
        finalizeCluster(std::move(cluster));
    }

    const size_t totalClusterCount = includedClusterCount + excludedClusterCount;
    Logger::log(LogLevel::INFO, "(", sampleName, ") Processed ", totalClusterCount,
                " clusters. Included ", includedClusterCount, " clusters, excluded ",
                excludedClusterCount, " clusters");

    return std::move(finishedClusters);
}

void InteractionClusterGenerator::finalizeCluster(InteractionCluster cluster) noexcept {
    logClusteringStatus();

    if (cluster.count < int(minReadCount)) {
        excludedClusterCount++;
        return;
    }

    includedClusterCount++;

    finishedClusters.emplace_back(std::move(cluster));
}

bool InteractionClusterGenerator::clusterIsBeforeOpenCluster(
    const InteractionCluster& cluster, const InteractionCluster& openCluster) const noexcept {
    return (cluster.segments.second.referenceIDIndex <
            openCluster.segments.second.referenceIDIndex) ||
           (cluster.segments.second.end < openCluster.segments.second.start);
}

void InteractionClusterGenerator::logClusteringStatus() const noexcept {
    const size_t totalClusterCount = includedClusterCount + excludedClusterCount;
    if (totalClusterCount % 100000 == 0 && finishedClusters.size() != 0) {
        Logger::log(LogLevel::INFO, "(", sampleName, ") Processed ", totalClusterCount,
                    " clusters. Included ", includedClusterCount, " clusters, excluded ",
                    excludedClusterCount, " clusters");
    }
}

}  // namespace analyze
}  // namespace pipelines
