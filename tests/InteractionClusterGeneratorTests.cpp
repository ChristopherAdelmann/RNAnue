#include <gtest/gtest.h>

#include <vector>

#include "DataTypes.hpp"
#include "InteractionClusterGenerator.hpp"
#include "gtest/gtest.h"

using namespace pipelines::analyze;

struct TestParam {
    const dtp::SplitRecords splitRecords;
    const std::vector<std::string> expectedSortOrderByID;
};

class InteractionClusterGeneratorTest : public testing::TestWithParam<TestParam> {
   protected:
    InteractionClusterGeneratorTest() : interactionClusterGenerator() {}

    InteractionClusterGenerator interactionClusterGenerator;
};

// TEST_P(InteractionClusterGeneratorTest, SplitRecordsAreSortedCorrectly) {
//     const auto& param = GetParam();
//     const auto& splitRecords = param.splitRecords;
//     const auto& expectedSortOrderByID = param.expectedSortOrderByID;

//     const auto& splitRecordsSorted = interactionClusterGenerator.sortSplitRecords(splitRecords);

//     for (size_t i = 0; i < splitRecordsSorted.size(); ++i) {
//         EXPECT_EQ(splitRecordsSorted[i].id, expectedSortOrderByID[i]);
//     }
// }
