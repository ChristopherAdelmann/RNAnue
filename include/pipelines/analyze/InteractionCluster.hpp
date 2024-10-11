#pragma once

// Standard
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Classes
#include "Segment.hpp"
#include "Utility.hpp"

namespace pipelines {
namespace analyze {

struct InteractionCluster {
    InteractionCluster(std::pair<Segment, Segment> segments,
                       std::vector<double> complementarityScores,
                       std::vector<double> hybridizationEnergies, int count = 1)
        : segments(std::move(segments)),
          complementarityScores(std::move(complementarityScores)),
          hybridizationEnergies(std::move(hybridizationEnergies)),
          count(count) {}

    std::pair<Segment, Segment> segments;
    std::vector<double> complementarityScores;
    std::vector<double> hybridizationEnergies;
    int count;
    std::optional<std::pair<std::string, std::string>> transcriptIDs = std::nullopt;
    std::optional<double> pValue = std::nullopt;
    std::optional<double> pAdj = std::nullopt;

    static std::optional<InteractionCluster> fromSegments(const Segment &segment1,
                                                          const Segment &segment2);

    bool operator<(const InteractionCluster &a) const;
    bool operator>(const InteractionCluster &a) const;
    bool operator==(const InteractionCluster &a) const;

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

inline std::ostream &operator<<(std::ostream &os, const InteractionCluster &interactionCluster) {
    return os << "InteractionCluster:\n"
              << "First segment: " << interactionCluster.segments.first
              << "\nSecond segment id: " << interactionCluster.segments.second << "\n"
              << "Complementarity score count: " << interactionCluster.complementarityScores.size()
              << "\n"
              << "Hybridization energie count: " << interactionCluster.hybridizationEnergies.size()
              << "\n"
              << "Count: " << interactionCluster.count << "\n";
};

}  // namespace analyze
}  // namespace pipelines
