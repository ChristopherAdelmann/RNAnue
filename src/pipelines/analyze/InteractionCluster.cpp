#include "InteractionCluster.hpp"

#include <ostream>

#include "Utility.hpp"

namespace pipelines {
namespace analyze {

InteractionCluster::InteractionCluster(std::pair<Segment, Segment> segments,
                                       const std::vector<double> &complementarityScores,
                                       const std::vector<double> &hybridizationEnergies)
    : segments(segments),
      complementarityScores(complementarityScores),
      hybridizationEnergies(hybridizationEnergies) {}

std::optional<std::pair<Segment, Segment>> InteractionCluster::getSortedElements(
    const Segment &segment1, const Segment &segment2) {
    if (segment1.recordID != segment2.recordID) [[unlikely]] {
        Logger::log(LogLevel::WARNING, "Record IDs do not match: ", segment1.recordID, " vs. ",
                    segment2.recordID, ". Make sure reads are sorted by name.");
        return std::nullopt;
    }

    return segment1.referenceIDIndex < segment2.referenceIDIndex ||
                   (segment1.referenceIDIndex == segment2.referenceIDIndex &&
                    segment1.start < segment2.start)
               ? std::make_pair(segment1, segment2)
               : std::make_pair(segment2, segment1);
}

std::optional<InteractionCluster> InteractionCluster::fromSegments(const Segment &segment1,
                                                                   const Segment &segment2) {
    const auto sortedElements = getSortedElements(segment1, segment2);
    if (!sortedElements.has_value()) [[unlikely]] {
        return std::nullopt;
    }

    assert(sortedElements.value().first.maxComplementarityScore ==
           sortedElements.value().second.maxComplementarityScore);
    assert(sortedElements.value().first.minHybridizationEnergy ==
           sortedElements.value().second.minHybridizationEnergy);

    return InteractionCluster(sortedElements.value(),
                              {sortedElements->first.maxComplementarityScore},
                              {sortedElements->first.minHybridizationEnergy}, 1);
}

bool InteractionCluster::operator<(const InteractionCluster &a) const {
    return std::tie(segments.first.start, segments.second.start) <
           std::tie(a.segments.first.start, a.segments.second.start);
}

bool InteractionCluster::operator>(const InteractionCluster &a) const {
    return std::tie(segments.first.start, segments.second.start) >
           std::tie(a.segments.first.start, a.segments.second.start);
}

bool InteractionCluster::operator==(const InteractionCluster &a) const {
    auto compareVectors = [](const std::vector<double> &v1, const std::vector<double> &v2) {
        return std::all_of(v1.begin(), v1.end(), [&v2](double val) {
            return std::any_of(v2.begin(), v2.end(),
                               [val](double a_val) { return helper::isEqual(val, a_val, 1e-6); });
        });
    };

    return compareVectors(complementarityScores, a.complementarityScores) &&
           compareVectors(hybridizationEnergies, a.hybridizationEnergies) &&
           segments.first == a.segments.first && segments.second == a.segments.second &&
           count == a.count;
}

bool InteractionCluster::overlaps(const InteractionCluster &other, const int graceDistance) const {
    bool isSameReferenceAndStrand =
        segments.first.referenceIDIndex == other.segments.first.referenceIDIndex &&
        segments.second.referenceIDIndex == other.segments.second.referenceIDIndex &&
        segments.first.strand == other.segments.first.strand &&
        segments.second.strand == other.segments.second.strand;

    if (!isSameReferenceAndStrand) {
        return false;
    }

    const bool firstOverlaps = segments.first.end + graceDistance >= other.segments.first.start &&
                               segments.first.start <= other.segments.first.end + graceDistance;
    const bool secondOverlaps =
        segments.second.end + graceDistance >= other.segments.second.start &&
        segments.second.start <= other.segments.second.end + graceDistance;

    return firstOverlaps && secondOverlaps;
}

void InteractionCluster::merge(const InteractionCluster &other) {
    assert(segments.first.referenceIDIndex == other.segments.first.referenceIDIndex);
    assert(segments.second.referenceIDIndex == other.segments.second.referenceIDIndex);
    assert(segments.first.strand == other.segments.first.strand);
    assert(segments.second.strand == other.segments.second.strand);

    segments.first.merge(other.segments.first);
    segments.second.merge(other.segments.second);

    complementarityScores.reserve(complementarityScores.size() +
                                  other.complementarityScores.size());
    hybridizationEnergies.reserve(hybridizationEnergies.size() +
                                  other.hybridizationEnergies.size());

    complementarityScores.insert(complementarityScores.end(),
                                 std::make_move_iterator(other.complementarityScores.begin()),
                                 std::make_move_iterator(other.complementarityScores.end()));
    hybridizationEnergies.insert(hybridizationEnergies.end(),
                                 std::make_move_iterator(other.hybridizationEnergies.begin()),
                                 std::make_move_iterator(other.hybridizationEnergies.end()));

    count += other.count;
}

double InteractionCluster::complementarityStatistics() const {
    assert(!complementarityScores.empty());
    assert(segments.first.maxComplementarityScore == segments.second.maxComplementarityScore);
    return helper::calculateMedian(complementarityScores) * segments.first.maxComplementarityScore;
}

double InteractionCluster::hybridizationEnergyStatistics() const {
    assert(!hybridizationEnergies.empty());
    assert(segments.first.minHybridizationEnergy == segments.second.minHybridizationEnergy);
    return std::sqrt(helper::calculateMedian(hybridizationEnergies) *
                     segments.first.minHybridizationEnergy);
}

}  // namespace analyze
}  // namespace pipelines
