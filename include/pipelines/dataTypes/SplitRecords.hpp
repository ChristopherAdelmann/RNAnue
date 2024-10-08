#pragma once

// Standard
#include <vector>

// Classes
#include "SamRecord.hpp"

namespace pipelines {
namespace dataTypes {

// Represents a single chimeric read which is split up into its blocks
struct SplitRecords {
    SplitRecords(const std::vector<SamRecord>& records) : records(records) {}

    std::vector<SamRecord> records;

    inline bool operator<(const SplitRecords& rhs) const {
        const auto& lhs_ref_id = records.back().reference_id();
        const auto& rhs_ref_id = rhs.records.back().reference_id();

        if (lhs_ref_id != rhs_ref_id) {
            return lhs_ref_id < rhs_ref_id;
        }

        const auto lhsEnd = recordEndPosition(records.back());
        const auto rhsEnd = recordEndPosition(rhs.records.back());

        return lhsEnd < rhsEnd;
    }

    inline bool operator>(const SplitRecords& rhs) const { return rhs < *this; }
};

// Split records are sorted by the end position of the last split record

}  // namespace dataTypes
}  // namespace pipelines
