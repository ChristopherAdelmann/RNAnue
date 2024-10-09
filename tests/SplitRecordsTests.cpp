#include <gtest/gtest.h>

// Standard
#include <algorithm>
#include <functional>
#include <ostream>
#include <string>
#include <vector>

// Classes
#include "ParseSamRecords.hpp"
#include "SplitRecords.hpp"

using namespace dataTypes;

struct SplitRecordsTestParams {
    const std::vector<SplitRecords> splitRecords;
    const std::vector<std::string> expectedBackRecordIDOrder;
};

class SplitRecordsTests : public ::testing::TestWithParam<SplitRecordsTestParams> {};

TEST_P(SplitRecordsTests, IsSortedFromBackToFront) {
    const SplitRecordsTestParams& param = GetParam();
    std::vector<SplitRecords> splitRecordGroups = param.splitRecords;
    const std::vector<std::string>& expectedBackRecordIDOrder = param.expectedBackRecordIDOrder;

    EXPECT_EQ(expectedBackRecordIDOrder.size(), splitRecordGroups.size());

    sort(splitRecordGroups.begin(), splitRecordGroups.end(), std::greater());

    std::vector<std::string> backRecordIDOrder;
    for (const auto& splitRecords : splitRecordGroups) {
        backRecordIDOrder.push_back(splitRecords.back().id());
    }

    EXPECT_EQ(expectedBackRecordIDOrder, backRecordIDOrder);
}

auto splitRecords1 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.231	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	40	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

auto splitRecords2 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.233	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.234	0	chromosome1	50	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

auto splitRecords3 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.235	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.236	0	chromosome1	50	20	7M	*	0	0	ATCGCGT	@@@@@@@	AS:i:0	XS:i:0
)";

auto splitRecords4 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.230	0	chromosome2	40	20	7M	*	0	0	ATCGCGT	@@@@@@@	AS:i:0	XS:i:0
SRR18331301.229	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

void PrintTo(const SplitRecordsTestParams& param, std::ostream* os) {
    *os << "SplitRecordsTestParams{" << "\n" << "splitRecords back IDs: \n";

    for (const auto& splitRecords : param.splitRecords) {
        *os << "\t" << splitRecords.back().id()
            << " reference id: " << splitRecords.back().reference_id().value_or(-1) << "\n";
    }

    *os << "expectedBackRecordIDOrder: \n";

    for (const auto& recordID : param.expectedBackRecordIDOrder) {
        *os << recordID << "\n";
    }

    *os << "}" << "\n";
};

const SplitRecordsTestParams sameChromosomeTestParams =
    SplitRecordsTestParams{{parseSamRecords(splitRecords2), parseSamRecords(splitRecords3),
                            parseSamRecords(splitRecords1)},
                           {"SRR18331301.236", "SRR18331301.234", "SRR18331301.232"}};

const SplitRecordsTestParams differentChromosomeTestParams = SplitRecordsTestParams{
    {parseSamRecords(splitRecords2), parseSamRecords(splitRecords4), parseSamRecords(splitRecords3),
     parseSamRecords(splitRecords1)},
    {"SRR18331301.230", "SRR18331301.236", "SRR18331301.234", "SRR18331301.232"}};

INSTANTIATE_TEST_SUITE_P(Default, SplitRecordsTests,
                         testing::Values(sameChromosomeTestParams, differentChromosomeTestParams));
