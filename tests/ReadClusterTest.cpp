#include <gtest/gtest.h>

#include "Cluster.hpp"

class ReadClusterStatisticsTest : public testing::Test {
   protected:
    ReadClusterStatisticsTest()
        : cluster(ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                            Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18})
                      .value()) {
        cluster.complementarityScores = {0.1, 0.2, 0.3, 0.4, 0.5};
        cluster.hybridizationEnergies = {-18.0, -20.0, -22.0, -24.0, -26.0};
    };

    ReadCluster cluster;
};

TEST_F(ReadClusterStatisticsTest, ComplementarityStatistics) {
    EXPECT_DOUBLE_EQ(cluster.complementarityStatistics(), 0.15);
}

TEST_F(ReadClusterStatisticsTest, HybridizationEnergyStatistics) {
    EXPECT_NEAR(cluster.hybridizationEnergyStatistics(), 19.8997487421, 1e-3);
}

// ReadCluster overlap tests
TEST(ReadClusterTest, OverlapsExact) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record1", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record2", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();

    EXPECT_TRUE(cluster1.overlaps(cluster2, 0));
}

TEST(ReadClusterTest, OverlapsGrace) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record1", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::FORWARD, 11, 20, 0.5, -18},
                                  Segment{"record2", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();

    EXPECT_TRUE(cluster1.overlaps(cluster2, 1));
}

TEST(ReadClusterTest, NoOverlapsFirstNotOverlapping) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record1", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::FORWARD, 11, 20, 0.5, -18},
                                  Segment{"record2", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();

    EXPECT_FALSE(cluster1.overlaps(cluster2, 0));
}

TEST(ReadClusterTest, NoOverlapsSecondNotOverlapping) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record1", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record2", 1, Strand::FORWARD, 11, 20, 0.5, -18})
            .value();

    EXPECT_FALSE(cluster1.overlaps(cluster2, 0));
}

TEST(ReadClusterTest, NoOverlapsDifferentReferenceID) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record1", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 2, Strand::FORWARD, 0, 10, 0.5, -18},
                                  Segment{"record2", 1, Strand::FORWARD, 0, 10, 0.5, -18})
            .value();

    EXPECT_FALSE(cluster1.overlaps(cluster2, 0));
}

// ReadCluster merge tests
TEST(ReadClusterTest, BasicMerge) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 5, 15, 0.5, -18},
                                  Segment{"record1", 0, Strand::FORWARD, 20, 25, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::FORWARD, 0, 10, 0.8, -15},
                                  Segment{"record2", 0, Strand::FORWARD, 20, 30, 0.8, -15})
            .value();

    cluster1.merge(cluster2);

    EXPECT_EQ(cluster1.count, 2);
    EXPECT_EQ(cluster1.complementarityScores.size(), 2);
    EXPECT_EQ(cluster1.hybridizationEnergies.size(), 2);

    EXPECT_EQ(cluster1.segments.first.start, 0);
    EXPECT_EQ(cluster1.segments.first.end, 15);

    EXPECT_EQ(cluster1.segments.second.start, 20);
    EXPECT_EQ(cluster1.segments.second.end, 30);

    EXPECT_EQ(cluster1.segments.first.referenceIDIndex, 0);
    EXPECT_EQ(cluster1.segments.second.referenceIDIndex, 0);

    EXPECT_EQ(cluster1.segments.first.hybridizationEnergy, -15);
    EXPECT_EQ(cluster1.segments.second.hybridizationEnergy, -15);

    EXPECT_EQ(cluster1.segments.first.complementarityScore, 0.8);
    EXPECT_EQ(cluster1.segments.second.complementarityScore, 0.8);
}

TEST(ReadClusterTest, MergeWithDifferentReferenceIDIndex) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 1, Strand::FORWARD, 5, 15, 0.5, -18},
                                  Segment{"record1", 0, Strand::FORWARD, 20, 25, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::FORWARD, 20, 30, 0.8, -15},
                                  Segment{"record2", 1, Strand::FORWARD, 0, 10, 0.8, -15})
            .value();

    cluster1.merge(cluster2);

    EXPECT_EQ(cluster1.count, 2);
    EXPECT_EQ(cluster1.complementarityScores.size(), 2);
    EXPECT_EQ(cluster1.hybridizationEnergies.size(), 2);

    EXPECT_EQ(cluster1.segments.first.start, 20);
    EXPECT_EQ(cluster1.segments.first.end, 30);

    EXPECT_EQ(cluster1.segments.second.start, 0);
    EXPECT_EQ(cluster1.segments.second.end, 15);

    EXPECT_EQ(cluster1.segments.first.referenceIDIndex, 0);
    EXPECT_EQ(cluster1.segments.second.referenceIDIndex, 1);

    EXPECT_EQ(cluster1.segments.first.hybridizationEnergy, -15);
    EXPECT_EQ(cluster1.segments.second.hybridizationEnergy, -15);

    EXPECT_EQ(cluster1.segments.first.complementarityScore, 0.8);
    EXPECT_EQ(cluster1.segments.second.complementarityScore, 0.8);

    EXPECT_EQ(cluster1.segments.first.strand, Strand::FORWARD);
    EXPECT_EQ(cluster1.segments.second.strand, Strand::FORWARD);
}

TEST(ReadClusterTest, MergeWithDifferentStrands) {
    ReadCluster cluster1 =
        ReadCluster::fromSegments(Segment{"record1", 0, Strand::FORWARD, 5, 15, 0.5, -18},
                                  Segment{"record1", 0, Strand::FORWARD, 20, 25, 0.5, -18})
            .value();
    ReadCluster cluster2 =
        ReadCluster::fromSegments(Segment{"record2", 0, Strand::REVERSE, 0, 10, 0.8, -15},
                                  Segment{"record2", 0, Strand::REVERSE, 20, 30, 0.8, -15})
            .value();

    EXPECT_DEATH(cluster1.merge(cluster2), ".*");
}
