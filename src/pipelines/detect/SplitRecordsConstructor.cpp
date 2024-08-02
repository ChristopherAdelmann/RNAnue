#include "SplitRecordsConstructor.hpp"

#include <vector>

#include "DataTypes.hpp"

dtp::SplitRecordsVariantGroups SplitRecordsConstructor::getSplitRecords(
    const std::vector<dtp::SamRecord> &recordGroup) {
    std::unordered_map<size_t, std::vector<dtp::SamRecord>> recordHitsGroup{};
}
