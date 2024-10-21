#include "InteractionCluster.hpp"

#include "Utility.hpp"

namespace pipelines::analyze {

auto InteractionCluster::getSortedElements(const Segment &segment1, const Segment &segment2)
    -> std::optional<std::pair<Segment, Segment>> {
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

auto InteractionCluster::fromSegments(const Segment &segment1, const Segment &segment2)
    -> std::optional<InteractionCluster> {
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

/**
 * @brief Less-than comparison operator for InteractionCluster.
 *
 * This operator compares two InteractionCluster objects based on the end position
 * of the second segment and the referenceIndexID. It returns true if the referenceIndexID is less
 * in the current object or the end position of the segment in the current object are
 * lexicographically less than those in the provided object.
 *
 * @param a The InteractionCluster object to compare with.
 * @return true if the current object is less than the provided object, false otherwise.
 */
auto InteractionCluster::operator<(const InteractionCluster &other) const -> bool {
    if (segments.second.referenceIDIndex == other.segments.second.referenceIDIndex) {
        return segments.second.end < other.segments.second.end;
    } else {
        return segments.second.referenceIDIndex < other.segments.second.referenceIDIndex;
    }
}

bool InteractionCluster::operator>(const InteractionCluster &other) const {
    if (segments.second.referenceIDIndex == other.segments.second.referenceIDIndex) {
        return segments.second.end > other.segments.second.end;
    } else {
        return segments.second.referenceIDIndex > other.segments.second.referenceIDIndex;
    }
}

auto InteractionCluster::operator==(const InteractionCluster &other) const -> bool {
    constexpr double EPSILON = 1e-6;
    auto compareVectors = [](const std::vector<double> &vec1, const std::vector<double> &vec2) {
        return std::all_of(vec1.begin(), vec1.end(), [&vec2](double val) {
            return std::any_of(vec2.begin(), vec2.end(), [val](double a_val) {
                return helper::isEqual(val, a_val, EPSILON);
            });
        });
    };

    return compareVectors(complementarityScores, other.complementarityScores) &&
           compareVectors(hybridizationEnergies, other.hybridizationEnergies) &&
           segments.first == other.segments.first && segments.second == other.segments.second &&
           count == other.count;
}

auto InteractionCluster::overlaps(const InteractionCluster &other,
                                  const int graceDistance) const -> bool {
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

auto InteractionCluster::complementarityStatistics() const -> double {
    assert(!complementarityScores.empty());
    assert(segments.first.maxComplementarityScore == segments.second.maxComplementarityScore);
    return helper::calculateMedian(complementarityScores) * segments.first.maxComplementarityScore;
}

auto InteractionCluster::hybridizationEnergyStatistics() const -> double {
    assert(!hybridizationEnergies.empty());
    assert(segments.first.minHybridizationEnergy == segments.second.minHybridizationEnergy);
    return std::sqrt(helper::calculateMedian(hybridizationEnergies) *
                     segments.first.minHybridizationEnergy);
}

}  // namespace pipelines::analyze
