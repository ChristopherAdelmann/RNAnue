#pragma once

// Standard
#include <cstdint>
#include <deque>
#include <optional>
#include <string>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SamRecord.hpp"

using seqan3::operator""_tag;
using namespace dataTypes;

namespace pipelines::analyze {

struct Segment {
    std::string recordID;
    int32_t referenceIDIndex;
    dataTypes::Strand strand;
    int32_t start;
    int32_t end;
    double maxComplementarityScore;
    double minHybridizationEnergy;

    static auto fromSamRecord(const SamRecord &record) -> std::optional<Segment>;

    [[nodiscard]] auto toGenomicRegion(const std::deque<std::string> &referenceIDs) const
        -> dataTypes::GenomicRegion;

    [[nodiscard]] auto toFeature(
        const std::deque<std::string> &referenceIDs, const std::string &featureID,
        const std::string &featureType = "transcript") const -> dataTypes::GenomicFeature;

    void merge(const Segment &other);

    auto operator==(const Segment &other) const -> bool;
};

inline auto operator<<(std::ostream &ostream, const Segment &segment) -> std::ostream & {
    return ostream << "Segment{recordID: " << segment.recordID
                   << ", referenceIDIndex: " << segment.referenceIDIndex
                   << ", strand: " << segment.strand << ", start: " << segment.start
                   << ", end: " << segment.end
                   << ", maxComplementarityScore: " << segment.maxComplementarityScore
                   << ", minHybridizationEnergy: " << segment.minHybridizationEnergy << "}";
};

}  // namespace pipelines::analyze
