#pragma once

// Standard
#include <cstdint>
#include <deque>
#include <optional>
#include <string>

// Classes
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"

namespace pipelines {
namespace analyze {

using seqan3::operator""_tag;

struct Segment {
    std::string recordID;
    int32_t referenceIDIndex;
    dtp::Strand strand;
    int32_t start;
    int32_t end;
    double maxComplementarityScore;
    double minHybridizationEnergy;

    static std::optional<Segment> fromSamRecord(const dtp::SamRecord &record);

    dtp::GenomicRegion toGenomicRegion(const std::deque<std::string> &referenceIDs) const;

    dtp::Feature toFeature(const std::deque<std::string> &referenceIDs,
                           const std::string &featureID,
                           const std::string &featureType = "transcript") const;

    void merge(const Segment &other);
};

}  // namespace analyze
}  // namespace pipelines