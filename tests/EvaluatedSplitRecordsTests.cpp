// #include <gtest/gtest.h>

// #include <seqan3/io/sam_file/input.hpp>
// #include <seqan3/utility/type_list/type_list.hpp>
// #include <sstream>

// #include "EvaluatedSplitRecords.hpp"
// #include "FeatureAnnotator.hpp"

// using namespace seqan3::literals;
// struct IsSplicedTestParam {
//     double minComplementarity = 0.5;
//     double minComplementarityFraction = 0.5;
//     double mfeThreshold = -10;
//     int spliceFilterTolerance;
//     Annotation::Orientation annotationOrientation;

//     dtp::SplitRecords splitRecords;
//     bool isSpliced;
// };

// class EvaluatedSplitRecordsTests : public testing::TestWithParam<IsSplicedTestParam> {
//    protected:
//     EvaluatedSplitRecordsTests() = default;

//     const std::deque<std::string> referenceIDs = {"chromosome1"};

//     const dtp::FeatureMap featureMap = {{"chromosome1",
//                                          {{.referenceID = "chromosome1",
//                                            .type = "exon",
//                                            .startPosition = 1,
//                                            .endPosition = 10,
//                                            .strand = dtp::Strand::FORWARD,
//                                            .id = "exon1",
//                                            .groupID = "gene1"},
//                                           {.referenceID = "chromosome1",
//                                            .type = "exon",
//                                            .startPosition = 20,
//                                            .endPosition = 30,
//                                            .strand = dtp::Strand::FORWARD,
//                                            .id = "exon2",
//                                            .groupID = "gene1"},
//                                           {.referenceID = "chromosome1",
//                                            .type = "exon",
//                                            .startPosition = 40,
//                                            .endPosition = 50,
//                                            .strand = dtp::Strand::FORWARD,
//                                            .id = "exon3",
//                                            .groupID = "gene1"}}}};

//     Annotation::FeatureAnnotator featureAnnotator;
// };

// void PrintTo(const IsSplicedTestParam& param, std::ostream* os) {
//     *os << "minComplementarity: " << param.minComplementarity
//         << ", minComplementarityFraction: " << param.minComplementarityFraction
//         << ", mfeThreshold: " << param.mfeThreshold
//         << ", spliceFilterTolerance: " << param.spliceFilterTolerance
//         << ", annotationOrientation: " << param.annotationOrientation;
// };

// TEST_P(EvaluatedSplitRecordsTests, IsSplicedSplitRecord) {
//     const IsSplicedTestParam& param = GetParam();

//     EvaluatedSplitRecords::SplicingParameters splicingParameters{
//         {param.minComplementarity, param.minComplementarityFraction, param.mfeThreshold},
//         param.annotationOrientation,
//         param.spliceFilterTolerance,
//         referenceIDs};

//     const auto isSpliced = EvaluatedSplitRecords::isSplicedSplitRecord(
//         param.splitRecords, splicingParameters, featureAnnotator);

//     EXPECT_EQ(isSpliced, param.isSpliced);
// };

// std::vector<dtp::SamRecord> parseSamRecords(const char* samFileRaw) {
//     std::istringstream samStream(samFileRaw);

//     using sam_file_input_t =
//         seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, dtp::sam_field_ids,
//                                seqan3::type_list<seqan3::format_sam>>;
//     sam_file_input_t samFile(std::istringstream{samFileRaw}, seqan3::format_sam{});

//     std::vector<dtp::SamRecord> samRecordsVector;

//     for (auto& record : samFile) {
//         samRecordsVector.push_back(record);
//     }
//     return samRecordsVector;
// }

// auto noSpliceRaw = R"(
// @HD     VN:1.6
// @SQ     SN:chromosome1  LN:2961149
// seq1	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@
// AS:i:0	XS:i:0 seq2	0	chromosome1	10	20	5M	*	0	0
// ATCGC	@@@@@	AS:i:0	XS:i:0
// )";

// auto spliceRaw = R"(
// @HD     VN:1.6
// @SQ     SN:chromosome1  LN:2961149
// seq1	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@
// AS:i:0	XS:i:0 seq2	0	chromosome1	20	20	5M	*	0	0
// ATCGC	@@@@@	AS:i:0	XS:i:0
// )";

// INSTANTIATE_TEST_SUITE_P(
//     Default, EvaluatedSplitRecordsTests,
//     testing::Values(IsSplicedTestParam{.spliceFilterTolerance = 5,
//                                        .annotationOrientation = Annotation::Orientation::BOTH,
//                                        .splitRecords = parseSamRecords(noSpliceRaw),
//                                        .isSpliced = false},
//                     IsSplicedTestParam{.spliceFilterTolerance = 5,
//                                        .annotationOrientation = Annotation::Orientation::BOTH,
//                                        .splitRecords = parseSamRecords(spliceRaw),
//                                        .isSpliced = true}));
