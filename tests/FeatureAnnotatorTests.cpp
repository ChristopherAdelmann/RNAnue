#include <gtest/gtest.h>

#include "FeatureAnnotator.hpp"

using namespace Annotation;

struct TestParam {
    TestParam(dtp::GenomicRegion region, Annotation::Orientation orientation,
              std::vector<std::string> expectedFeatureIds)
        : region(region),
          expectedFeatureIds(std::move(expectedFeatureIds)),
          orientation(orientation) {}

    TestParam(dtp::GenomicRegion region, std::vector<std::string> expectedFeatureIds)
        : region(region), expectedFeatureIds(std::move(expectedFeatureIds)) {}

    dtp::GenomicRegion region;
    std::vector<std::string> expectedFeatureIds;
    Annotation::Orientation orientation = Annotation::Orientation::SAME;
};

// Tests for FeatureAnnotator::overlappingFeatures and FeatureAnnotator::overlappingFeatureIterator
class FeatureAnnotatorTest : public testing::TestWithParam<TestParam> {
   protected:
    FeatureAnnotatorTest() : annotator(featureMap) {}

    const dtp::FeatureMap featureMap = {
        {"chromosome1",
         {{"chromosome1", "transcript", 1, 10, dtp::Strand::FORWARD, "feature1"},
          {"chromosome1", "transcript", 20, 30, dtp::Strand::FORWARD, "feature2"},
          {"chromosome1", "transcript", 40, 50, dtp::Strand::FORWARD, "feature3"},
          {"chromosome1", "transcript", 5, 25, dtp::Strand::REVERSE, "feature4"}}},
        {"chromosome2",
         {{"chromosome2", "transcript", 1, 10, dtp::Strand::FORWARD, "feature3"},
          {"chromosome2", "transcript", 20, 30, dtp::Strand::FORWARD, "feature4"}}},
    };

    FeatureAnnotator annotator;
};

void PrintTo(const TestParam& param, std::ostream* os) {
    *os << "GenomicRegion{referenceID: " << param.region.referenceID
        << ", startPosition: " << param.region.startPosition
        << ", endPosition: " << param.region.endPosition << ", strand: "
        << (param.region.strand.has_value() ? std::to_string(param.region.strand.value())
                                            : "nullopt")
        << "}";
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatures) {
    const auto& param = GetParam();
    const auto features = annotator.overlappingFeatures(param.region, param.orientation);

    ASSERT_EQ(features.size(), param.expectedFeatureIds.size());
    for (size_t i = 0; i < features.size(); ++i) {
        EXPECT_EQ(features[i].id, param.expectedFeatureIds[i]);
    }
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatureIterator) {
    const auto& param = GetParam();
    const auto results = annotator.overlappingFeatureIterator(param.region, param.orientation);

    size_t i = 0;
    for (const auto& feature : results) {
        EXPECT_EQ(feature.id, param.expectedFeatureIds[i]);
        ++i;
    }
    EXPECT_EQ(i, param.expectedFeatureIds.size());
}

INSTANTIATE_TEST_SUITE_P(
    Default, FeatureAnnotatorTest,
    testing::Values(TestParam{{"chromosome1", 9, 15, dtp::Strand::FORWARD}, {"feature1"}},
                    TestParam{{"chromosome1", 5, 15, dtp::Strand::REVERSE}, {"feature4"}},
                    TestParam{{"chromosome1", 5, 15, std::nullopt}, {"feature1", "feature4"}},
                    TestParam{{"chromosome1", 31, 39, std::nullopt}, {}},
                    TestParam{{"chromosome2", 5, 15, std::nullopt}, {"feature3"}},
                    TestParam{{"chromosome3", 5, 25, std::nullopt}, {}},
                    TestParam{{"chromosome1", 5, 15, dtp::Strand::FORWARD},
                              Annotation::Orientation::OPPOSITE,
                              {"feature4"}}));

// Tests for FeatureAnnotator::getBestOverlappingFeature
class BestFeatureAnnotatorTest : public testing::TestWithParam<TestParam> {
   protected:
    BestFeatureAnnotatorTest() : annotator(featureMap) {}

    const dtp::FeatureMap featureMap = {
        {"chromosome1",
         {{"chromosome1", "transcript", 1, 10, dtp::Strand::FORWARD, "feature1"},
          {"chromosome1", "transcript", 20, 30, dtp::Strand::FORWARD, "feature2"},
          {"chromosome1", "transcript", 40, 50, dtp::Strand::FORWARD, "feature3"},
          {"chromosome1", "transcript", 5, 25, dtp::Strand::REVERSE, "feature4"}}},
        {"chromosome2",
         {{"chromosome2", "transcript", 1, 10, dtp::Strand::FORWARD, "feature3"},
          {"chromosome2", "transcript", 20, 30, dtp::Strand::FORWARD, "feature4"}}},
    };

    FeatureAnnotator annotator;
};

TEST_P(BestFeatureAnnotatorTest, GetBestOverlappingFeature) {
    const auto& param = GetParam();
    const auto feature =
        annotator.getBestOverlappingFeature(param.region, Annotation::Orientation::BOTH);

    if (param.expectedFeatureIds.empty()) {
        EXPECT_FALSE(feature.has_value());
    } else {
        EXPECT_TRUE(feature.has_value());
        EXPECT_EQ(feature->id, param.expectedFeatureIds[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    Default, BestFeatureAnnotatorTest,
    testing::Values(TestParam{{"chromosome1", 1, 7, dtp::Strand::FORWARD}, {"feature1"}},
                    TestParam{{"chromosome1", 5, 15, dtp::Strand::REVERSE}, {"feature4"}},
                    TestParam{{"chromosome1", 5, 15, std::nullopt}, {"feature4"}},
                    TestParam{{"chromosome1", 31, 39, std::nullopt}, {}},
                    TestParam{{"chromosome2", 5, 15, std::nullopt}, {"feature3"}},
                    TestParam{{"chromosome3", 5, 25, std::nullopt}, {}}));

// Tests for FeatureAnnotator::insert and FeatureAnnotator::mergeInsert
class InsertFeatureAnnotatorTest : public testing::Test {
   protected:
    InsertFeatureAnnotatorTest() : annotator(featureMap) {}

    const dtp::FeatureMap featureMap = {
        {"chromosome1",
         {{"chromosome1", "transcript", 1, 10, dtp::Strand::FORWARD, "feature1"},
          {"chromosome1", "transcript", 20, 30, dtp::Strand::FORWARD, "feature2"}}},
    };

    FeatureAnnotator annotator;
};

TEST_F(InsertFeatureAnnotatorTest, Insert) {
    const dtp::GenomicRegion region{"chromosome1", 5, 20, dtp::Strand::FORWARD};
    const auto featureId = annotator.insert(region);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 2);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_NE(features[1].id, "feature2");
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsert) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, dtp::Strand::FORWARD};
    const auto result = annotator.mergeInsert(region, 0);

    ASSERT_EQ(annotator.featureCount(), 2);

    EXPECT_EQ(result.featureID, "feature1");
    EXPECT_TRUE(result.mergedFeatureIDs.empty());

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertTwoOverlapping) {
    const dtp::GenomicRegion region{"chromosome1", 5, 25, dtp::Strand::FORWARD};
    const auto result = annotator.mergeInsert(region, 0);

    ASSERT_EQ(annotator.featureCount(), 1);

    EXPECT_EQ(result.featureID, "feature1");
    ASSERT_EQ(result.mergedFeatureIDs.size(), 1);
    EXPECT_EQ(result.mergedFeatureIDs[0], "feature2");

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 30);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertWithReverseStrand) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, dtp::Strand::REVERSE};
    const auto result = annotator.mergeInsert(region, 0);

    ASSERT_EQ(annotator.featureCount(), 3);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 5);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertWithNoStrand) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, std::nullopt};
    ASSERT_DEATH(annotator.mergeInsert(region, 0), "Strand must be specified for insertion");
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertNotExistingReferenceID) {
    const dtp::GenomicRegion region{"chromosome2", 5, 25, dtp::FORWARD};
    const auto result = annotator.mergeInsert(region, 0);

    ASSERT_EQ(annotator.featureCount(), 3);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 5);
    EXPECT_EQ(features[0].endPosition, 25);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistance) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, dtp::FORWARD};
    const auto result = annotator.mergeInsert(region, 5);

    ASSERT_EQ(annotator.featureCount(), 1);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_EQ(result.mergedFeatureIDs, std::vector<std::string>{"feature2"});
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 30);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistanceNotSecondOverlapping) {
    const dtp::GenomicRegion region{"chromosome1", 1, 14, dtp::FORWARD};
    const auto result = annotator.mergeInsert(region, 5);

    ASSERT_EQ(annotator.featureCount(), 2);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 14);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistanceOneSpace) {
    const dtp::GenomicRegion region{"chromosome1", 11, 15, dtp::FORWARD};
    const auto result = annotator.mergeInsert(region, 0);

    ASSERT_EQ(annotator.featureCount(), 3);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 11);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistanceBluntEnds) {
    const dtp::GenomicRegion region{"chromosome1", 10, 15, dtp::FORWARD};
    const auto result = annotator.mergeInsert(region, 0);

    ASSERT_EQ(annotator.featureCount(), 2);

    const auto features = annotator.overlappingFeatures(region, Annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 15);
}