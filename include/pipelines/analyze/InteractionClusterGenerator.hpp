#pragma once

// Standard
#include <algorithm>
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
    InteractionClusterGenerator() = default;
    ~InteractionClusterGenerator() = default;

    using InteractionClusters = std::vector<InteractionCluster>;

    InteractionClusters mergeClusters(InteractionClusters& clusters,
                                      const int graceDistance) noexcept;

   private:
    InteractionClusters finishedClusters;
    std::forward_list<InteractionCluster> openClusterQueue;

    void finalizeCluster(InteractionCluster cluster) noexcept;
    bool clusterIsBeforeOpenCluster(const InteractionCluster& cluster,
                                    const InteractionCluster& openCluster) const noexcept;
    void logClusteringStatus() const noexcept;
};
}  // namespace analyze
}  // namespace pipelines
