#include <gtest/gtest.h>

#include "FeatureAnnotator.hpp"

using namespace Annotation;

struct TestParam {
    dtp::GenomicRegion region;
    std::vector<std::string> expectedFeatureIds;
};

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
    const auto features = annotator.overlappingFeatures(param.region);

    ASSERT_EQ(features.size(), param.expectedFeatureIds.size());
    for (size_t i = 0; i < features.size(); ++i) {
        EXPECT_EQ(features[i].id, param.expectedFeatureIds[i]);
    }
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatureIterator) {
    const auto& param = GetParam();
    const auto results = annotator.overlappingFeatureIterator(param.region);

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
                    TestParam{{"chromosome3", 5, 25, std::nullopt}, {}}));

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
    const auto feature = annotator.getBestOverlappingFeature(param.region);

    if (param.expectedFeatureIds.empty()) {
        EXPECT_FALSE(feature.has_value());
    } else {
        EXPECT_TRUE(feature.has_value());
        EXPECT_EQ(feature->id, param.expectedFeatureIds[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    Default, BestFeatureAnnotatorTest,
    testing::Values(TestParam{{"chromosome1", 9, 15, dtp::Strand::FORWARD}, {"feature1"}},
                    TestParam{{"chromosome1", 5, 15, dtp::Strand::REVERSE}, {"feature4"}},
                    TestParam{{"chromosome1", 5, 15, std::nullopt}, {"feature4"}},
                    TestParam{{"chromosome1", 31, 39, std::nullopt}, {}},
                    TestParam{{"chromosome2", 5, 15, std::nullopt}, {"feature3"}},
                    TestParam{{"chromosome3", 5, 25, std::nullopt}, {}}));

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

    const auto features = annotator.overlappingFeatures(region);
    ASSERT_EQ(features.size(), 2);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_NE(features[1].id, "feature2");
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsert) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, dtp::Strand::FORWARD};
    const auto featureId = annotator.mergeInsert(region);

    ASSERT_EQ(annotator.featureCount(), 2);

    EXPECT_EQ(featureId, "feature1");

    const auto features = annotator.overlappingFeatures(region);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertTwoOverlapping) {
    const dtp::GenomicRegion region{"chromosome1", 5, 25, dtp::Strand::FORWARD};
    const auto featureId = annotator.mergeInsert(region);

    ASSERT_EQ(annotator.featureCount(), 1);

    EXPECT_EQ(featureId, "feature1");

    const auto features = annotator.overlappingFeatures(region);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 30);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertWithReverseStrand) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, dtp::Strand::REVERSE};
    const auto featureId = annotator.mergeInsert(region);

    ASSERT_EQ(annotator.featureCount(), 3);

    const auto features = annotator.overlappingFeatures(region);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, featureId);
    EXPECT_EQ(features[0].startPosition, 5);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertWithNoStrand) {
    const dtp::GenomicRegion region{"chromosome1", 5, 15, std::nullopt};
    ASSERT_DEATH(annotator.mergeInsert(region), "Strand must be specified for insertion");
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertNotExistingReferenceID) {
    const dtp::GenomicRegion region{"chromosome2", 5, 25, dtp::FORWARD};
    const auto featureId = annotator.mergeInsert(region);

    ASSERT_EQ(annotator.featureCount(), 3);

    const auto features = annotator.overlappingFeatures(region);
    ASSERT_EQ(features.size(), 1);
    EXPECT_EQ(features[0].id, featureId);
    EXPECT_EQ(features[0].startPosition, 5);
    EXPECT_EQ(features[0].endPosition, 25);
}