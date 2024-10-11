#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <execution>
#include <forward_list>
#include <functional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

// Classes
#include "InteractionCluster.hpp"
#include "Logger.hpp"

namespace pipelines {
namespace analyze {

class InteractionClusterGenerator {
   public:
    InteractionClusterGenerator(size_t minReadCount, int graceDistance) noexcept
        : minReadCount(minReadCount), graceDistance(graceDistance) {}
    ~InteractionClusterGenerator() = default;

    using InteractionClusters = std::vector<InteractionCluster>;

    InteractionClusters mergeClusters(InteractionClusters& clusters) noexcept;

   private:
    InteractionClusters finishedClusters;
    std::forward_list<InteractionCluster> openClusterQueue;

    size_t minReadCount;
    int graceDistance;

    size_t includedClusterCount = 0;
    size_t excludedClusterCount = 0;

    void finalizeCluster(InteractionCluster cluster) noexcept;
    bool clusterIsBeforeOpenCluster(const InteractionCluster& cluster,
                                    const InteractionCluster& openCluster) const noexcept;
    void logClusteringStatus() const noexcept;
};
}  // namespace analyze
}  // namespace pipelines
