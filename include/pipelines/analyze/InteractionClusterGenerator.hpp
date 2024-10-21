#pragma once

// Standard
#include <cstddef>
#include <forward_list>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "InteractionCluster.hpp"

namespace pipelines::analyze {

class InteractionClusterGenerator {
   public:
    InteractionClusterGenerator(std::string sampleName, size_t minReadCount,
                                int graceDistance) noexcept
        : sampleName(std::move(sampleName)),
          minReadCount(minReadCount),
          graceDistance(graceDistance) {}

    using InteractionClusters = std::vector<InteractionCluster>;

    auto mergeClusters(InteractionClusters& clusters) noexcept -> InteractionClusters;

   private:
    InteractionClusters finishedClusters;
    std::forward_list<InteractionCluster> openClusterQueue;

    std::string sampleName;

    size_t minReadCount;
    int graceDistance;

    size_t includedClusterCount = 0;
    size_t excludedClusterCount = 0;

    void finalizeCluster(InteractionCluster cluster) noexcept;
    static auto clusterIsBeforeOpenCluster(const InteractionCluster& cluster,
                                           const InteractionCluster& openCluster) noexcept -> bool;
    void logClusteringStatus() const noexcept;
};
}  // namespace pipelines::analyze
