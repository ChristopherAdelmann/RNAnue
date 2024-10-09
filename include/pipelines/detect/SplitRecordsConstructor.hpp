#pragma once

// Standard
#include <cstddef>
#include <unordered_map>

// Classes
#include <vector>

#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsVariantGroups.hpp"

using namespace dataTypes;

class SplitRecordsConstructor {
   public:
    SplitRecordsConstructor() = delete;
    ~SplitRecordsConstructor() = delete;

    static SplitRecordsVariantGroups getSplitRecords(const std::vector<SamRecord>& recordGroup);

   private:
    static SplitRecords constructSplitRecordsForHitGroup(const std::vector<SamRecord>& hitGroup);

    static SplitRecords constructSplitRecords(const SamRecord& record);
};
