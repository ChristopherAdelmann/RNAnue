#pragma once

// Classes
#include <filesystem>
#include <vector>

#include "InteractionCluster.hpp"

namespace pipelines {
namespace analyze {

namespace fs = std::filesystem;

struct SplitReadParser {
    static std::vector<InteractionCluster> parse(const fs::path splitRecordsFilePath) {
        seqan3::sam_file_input splitsIn{splitRecordsFilePath, sam_field_ids{}};

        std::vector<InteractionCluster> clusters;

        for (auto &&records : splitsIn | seqan3::views::chunk(2)) {
            std::optional<Segment> segment1 = Segment::fromSamRecord(*records.begin());
            std::optional<Segment> segment2 = Segment::fromSamRecord(*(++records.begin()));

            if (!segment1 || !segment2) [[unlikely]] {
                continue;
            }

            auto cluster = InteractionCluster::fromSegments(*segment1, *segment2);

            if (!cluster) [[unlikely]] {
                continue;
            }

            clusters.push_back(*cluster);
        }

        return clusters;
    };
};
}  // namespace analyze
}  // namespace pipelines