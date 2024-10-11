#include "InteractionClusterGenerator.hpp"

#include <ranges>
#include <string>
#include <utility>

#include "Logger.hpp"

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
    InteractionClusterGenerator::InteractionClusters& clusters, const int graceDistance) noexcept {
    // Clusters should be sorted from back to front
    std::sort(clusters.begin(), clusters.end(), std::greater<>{});

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

    Logger::log(LogLevel::INFO, "Found ", std::to_string(finishedClusters.size()), " clusters");

    return std::move(finishedClusters);
}

void InteractionClusterGenerator::finalizeCluster(InteractionCluster cluster) noexcept {
    logClusteringStatus();

    finishedClusters.emplace_back(std::move(cluster));
}

bool InteractionClusterGenerator::clusterIsBeforeOpenCluster(
    const InteractionCluster& cluster, const InteractionCluster& openCluster) const noexcept {
    return (cluster.segments.second.referenceIDIndex <
            openCluster.segments.second.referenceIDIndex) ||
           (cluster.segments.second.end < openCluster.segments.second.start);
}

void InteractionClusterGenerator::logClusteringStatus() const noexcept {
    if (finishedClusters.size() % 50 == 0 && finishedClusters.size() != 0) {
        Logger::log(LogLevel::INFO, "Found ", std::to_string(finishedClusters.size()), " clusters");
    }
}

}  // namespace analyze
}  // namespace pipelines
