#pragma once

// Standard
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "Segment.hpp"

namespace pipelines::analyze {

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

    static auto fromSegments(const Segment &segment1,
                             const Segment &segment2) -> std::optional<InteractionCluster>;

    auto operator<(const InteractionCluster &other) const -> bool;
    auto operator>(const InteractionCluster &other) const -> bool;
    auto operator==(const InteractionCluster &other) const -> bool;

    [[nodiscard]] auto overlaps(const InteractionCluster &other, int graceDistance) const -> bool;

    void merge(const InteractionCluster &other);

    [[nodiscard]] auto complementarityStatistics() const -> double;

    [[nodiscard]] auto hybridizationEnergyStatistics() const -> double;

   private:
    InteractionCluster(std::pair<Segment, Segment> segments,
                       const std::vector<double> &complementarityScores,
                       const std::vector<double> &hybridizationEnergies);

    static auto getSortedElements(const Segment &segment1, const Segment &segment2)
        -> std::optional<std::pair<Segment, Segment>>;
};

inline auto operator<<(std::ostream &os,
                       const InteractionCluster &interactionCluster) -> std::ostream & {
    return os << "InteractionCluster:\n"
              << "First segment: " << interactionCluster.segments.first
              << "\nSecond segment id: " << interactionCluster.segments.second << "\n"
              << "Complementarity score count: " << interactionCluster.complementarityScores.size()
              << "\n"
              << "Hybridization energie count: " << interactionCluster.hybridizationEnergies.size()
              << "\n"
              << "Count: " << interactionCluster.count << "\n";
};

}  // namespace pipelines::analyze
