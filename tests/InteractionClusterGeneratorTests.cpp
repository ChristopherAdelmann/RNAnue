#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "DataTypes.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "Segment.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsParser.hpp"
#include "gtest/gtest.h"
using namespace pipelines::analyze;

// Note: This is the expected clustering result:
// Cluster 1: 1,3,6
// Cluster 2: 2,7
// Cluster 3: 4,5

using namespace pipelines::analyze;

std::string testInteractionClustersSamPath() {
    return (std::filesystem::path{__FILE__}.parent_path() /
            "test_data/interactionClusterRecords.sam")
        .string();
}

struct InteractionClusterGeneratorTestParam {
    const std::vector<InteractionCluster> expectedInteractionClusters;
};

class InteractionClusterGeneratorTests
    : public testing::TestWithParam<InteractionClusterGeneratorTestParam> {
   protected:
    InteractionClusterGeneratorTests() : interactionClusterGenerator(0, 1) {}

    InteractionClusterGenerator interactionClusterGenerator;
};

TEST_P(InteractionClusterGeneratorTests, SplitRecordsAreSortedCorrectly) {
    const auto& param = GetParam();
    const auto& expectedClusters = param.expectedInteractionClusters;

    std::vector<InteractionCluster> interactionClusters =
        SplitReadParser::parse(testInteractionClustersSamPath());

    EXPECT_EQ(interactionClusters.size(), 7u);

    std::vector<InteractionCluster> mergedInteractionClusters =
        interactionClusterGenerator.mergeClusters(interactionClusters);

    EXPECT_EQ(mergedInteractionClusters.size(), expectedClusters.size());

    for (const auto& mergedCluster : mergedInteractionClusters) {
        auto it = std::find_if(expectedClusters.begin(), expectedClusters.end(),
                               [&mergedCluster](const InteractionCluster& expectedCluster) {
                                   return mergedCluster == expectedCluster;
                               });

        EXPECT_NE(it, expectedClusters.end());
    }
}

const InteractionCluster cluster1 =
    InteractionCluster(std::pair{Segment{.recordID = "SRR18331301.3",
                                         .referenceIDIndex = 0,
                                         .strand = dtp::Strand::FORWARD,
                                         .start = 19,
                                         .end = 27,
                                         .maxComplementarityScore = 1,
                                         .minHybridizationEnergy = -1.7},
                                 Segment{.recordID = "SRR18331301.3",
                                         .referenceIDIndex = 1,
                                         .strand = dtp::Strand::FORWARD,
                                         .start = 49,
                                         .end = 64,
                                         .maxComplementarityScore = 1,
                                         .minHybridizationEnergy = -1.7}},
                       {1, 1, 1}, {-1.7, -1.7, -1.7}, 3);

const InteractionCluster cluster2 =
    InteractionCluster(std::pair{Segment{.recordID = "SRR18331301.2",
                                         .referenceIDIndex = 1,
                                         .strand = dtp::Strand::FORWARD,
                                         .start = 4,
                                         .end = 10,
                                         .maxComplementarityScore = 1,
                                         .minHybridizationEnergy = -1.7},
                                 Segment{.recordID = "SRR18331301.2",
                                         .referenceIDIndex = 1,
                                         .strand = dtp::Strand::FORWARD,
                                         .start = 51,
                                         .end = 57,
                                         .maxComplementarityScore = 1,
                                         .minHybridizationEnergy = -1.7}},
                       {1, 1}, {-1.7, -1.7}, 2);

const InteractionCluster cluster3 =
    InteractionCluster(std::pair{Segment{.recordID = "SRR18331301.4",
                                         .referenceIDIndex = 0,
                                         .strand = dtp::Strand::FORWARD,
                                         .start = 4,
                                         .end = 14,
                                         .maxComplementarityScore = 1,
                                         .minHybridizationEnergy = -1.7},
                                 Segment{.recordID = "SRR18331301.4",
                                         .referenceIDIndex = 0,
                                         .strand = dtp::Strand::FORWARD,
                                         .start = 39,
                                         .end = 46,
                                         .maxComplementarityScore = 1,
                                         .minHybridizationEnergy = -1.7}},
                       {1, 1}, {-1.7, -1.7}, 2);

INSTANTIATE_TEST_SUITE_P(Default, InteractionClusterGeneratorTests,
                         testing::Values(InteractionClusterGeneratorTestParam{
                             {cluster1, cluster2, cluster3}}));
