#pragma once

// Standard
#include <cstddef>
#include <unordered_map>

// Classes
#include <vector>

#include "DataTypes.hpp"

class SplitRecordsConstructor {
   public:
    SplitRecordsConstructor() = delete;
    ~SplitRecordsConstructor() = delete;

    static dtp::SplitRecordsVariantGroups getSplitRecords(
        const std::vector<dtp::SamRecord>& recordGroup);

   private:
    static dtp::SplitRecords constructSplitRecordsForHitGroup(
        const std::vector<dtp::SamRecord>& hitGroup);

    static dtp::SplitRecords constructSplitRecords(const dtp::SamRecord& record);
};
