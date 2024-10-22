#include "Segment.hpp"

// Standard
#include <optional>

// Include
#include "CustomSamTags.hpp"
#include "GenomicRegion.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {

auto Segment::fromSamRecord(const SamRecord &record) -> std::optional<Segment> {
    if (!record.reference_position().has_value() || !record.reference_id().has_value()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const dataTypes::Strand strand{isReverseStrand ? dataTypes::Strand::REVERSE
                                                   : dataTypes::Strand::FORWARD};

    const auto start = record.reference_position();
    const std::optional<int32_t> end = dataTypes::recordEndPosition(record);

    if (!start.has_value() || !end.has_value()) {
        return std::nullopt;
    }

    const double hybridizationEnergy = record.tags().get<"XE"_tag>();
    const double complementarityScore = record.tags().get<"XC"_tag>();

    return Segment{record.id(),
                   record.reference_id().value(),
                   strand,
                   start.value(),
                   end.value(),
                   complementarityScore,
                   hybridizationEnergy};
}

auto Segment::toGenomicRegion(const std::deque<std::string> &referenceIDs) const
    -> dataTypes::GenomicRegion {
    return dataTypes::GenomicRegion{referenceIDs[referenceIDIndex], start, end, strand};
}

auto Segment::toFeature(const std::deque<std::string> &referenceIDs, const std::string &featureID,
                        const std::string &featureType) const -> dataTypes::GenomicFeature {
    return dataTypes::GenomicFeature{referenceIDs[referenceIDIndex],
                                     featureType,
                                     start,
                                     end,
                                     strand,
                                     featureID,
                                     std::nullopt,
                                     std::nullopt};
}

void Segment::merge(const Segment &other) {
    start = std::min(start, other.start);
    end = std::max(end, other.end);
    maxComplementarityScore = std::max(maxComplementarityScore, other.maxComplementarityScore);
    minHybridizationEnergy = std::min(minHybridizationEnergy, other.minHybridizationEnergy);
}

auto Segment::operator==(const Segment &other) const -> bool {
    constexpr double EPSILON = 1e-6;
    return referenceIDIndex == other.referenceIDIndex && strand == other.strand &&
           start == other.start && end == other.end &&
           helper::isEqual(maxComplementarityScore, other.maxComplementarityScore, EPSILON) &&
           helper::isEqual(minHybridizationEnergy, other.minHybridizationEnergy, EPSILON);
}

}  // namespace pipelines::analyze
