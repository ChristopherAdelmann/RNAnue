#pragma once

// Standard
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Classes
#include "Segment.hpp"
#include "Utility.hpp"

namespace pipelines {
namespace analyze {

struct InteractionCluster {
    std::pair<Segment, Segment> segments;
    std::vector<double> complementarityScores;
    std::vector<double> hybridizationEnergies;
    std::optional<std::pair<std::string, std::string>> transcriptIDs = std::nullopt;
    std::optional<double> pValue = std::nullopt;
    std::optional<double> pAdj = std::nullopt;
    int count{1};

    static std::optional<InteractionCluster> fromSegments(const Segment &segment1,
                                                          const Segment &segment2);

    bool operator<(const InteractionCluster &a) const;

    bool overlaps(const InteractionCluster &other, const int graceDistance) const;

    void merge(const InteractionCluster &other);

    double complementarityStatistics() const;

    double hybridizationEnergyStatistics() const;

   private:
    InteractionCluster(std::pair<Segment, Segment> segments,
                       const std::vector<double> &complementarityScores,
                       const std::vector<double> &hybridizationEnergies);

    static std::optional<std::pair<Segment, Segment>> getSortedElements(const Segment &segment1,
                                                                        const Segment &segment2);
};

}  // namespace analyze
}  // namespace pipelines
