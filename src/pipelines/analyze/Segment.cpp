#include "Segment.hpp"

#include <optional>

namespace pipelines {
namespace analyze {

std::optional<Segment> Segment::fromSamRecord(const dtp::SamRecord &record) {
    if (!record.reference_position().has_value() || !record.reference_id().has_value()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const dtp::Strand strand{isReverseStrand ? dtp::Strand::REVERSE : dtp::Strand::FORWARD};

    const auto start = record.reference_position();
    const std::optional<int32_t> end = dtp::recordEndPosition(record);

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

dtp::GenomicRegion Segment::toGenomicRegion(const std::deque<std::string> &referenceIDs) const {
    return dtp::GenomicRegion{referenceIDs[referenceIDIndex], start, end, strand};
}

dtp::Feature Segment::toFeature(const std::deque<std::string> &referenceIDs,
                                const std::string &featureID,
                                const std::string &featureType) const {
    return dtp::Feature{referenceIDs[referenceIDIndex],
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

}  // namespace analyze
}  // namespace pipelines
