#pragma once

// Standard
#include <algorithm>
#include <deque>
#include <filesystem>
#include <fstream>
#include <future>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>

// Classes
#include "AsyncSplitRecordGroupBuffer.hpp"
#include "Constants.hpp"
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"
#include "DetectData.hpp"
#include "DetectParameters.hpp"
#include "DetectSample.hpp"
#include "FeatureAnnotator.hpp"
#include "Logger.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsEvaluator.hpp"
#include "Utility.hpp"

using namespace dataTypes;

namespace pipelines {
namespace detect {
namespace fs = std::filesystem;

class Detect {
   public:
    explicit Detect(DetectParameters params)
        : params(params),
          featureAnnotator(params.featuresInPath, params.featureTypes),
          splitRecordsEvaluator(SplitRecordsEvaluator(getSplitRecordsEvaluatorParameters(params))) {
          };
    ~Detect() = default;

    void process(const DetectData &data);

   private:
    using AsyncGroupBufferType = AsyncSplitRecordGroupBufferView<std::ranges::ref_view<
        seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, sam_field_ids>>>;

    struct ChunkedOutTmpDirs {
        fs::path outputTmpSplitsDir;
        fs::path outputTmpMultisplitsDir;
        fs::path outputTmpUnassignedContiguousDir;
    };

    using TranscriptCounts = std::unordered_map<std::string, size_t>;

    struct Result {
        size_t processedRecordsCount{0};
        TranscriptCounts transcriptCounts;
        size_t splitFragmentsCount{0};
        size_t singletonFragmentsCount{0};

        void operator+=(const Result &other) {
            processedRecordsCount += other.processedRecordsCount;
            splitFragmentsCount += other.splitFragmentsCount;
            singletonFragmentsCount += other.singletonFragmentsCount;

            for (const auto &[transcript, count] : other.transcriptCounts) {
                transcriptCounts[transcript] += count;
            }
        }
    };

    DetectParameters params;

    annotation::FeatureAnnotator featureAnnotator;
    SplitRecordsEvaluator splitRecordsEvaluator;

    const SplitRecordsEvaluationParameters::ParameterVariant getSplitRecordsEvaluatorParameters(
        const DetectParameters &params) const;

    const std::deque<std::string> &getReferenceIDs(const fs::path &mappingsInPath) const;

    void processSample(const DetectSample &sample) const;

    Result processRecordChunk(const ChunkedOutTmpDirs &outTmpDirs,
                              AsyncGroupBufferType &recordInputBuffer,
                              const std::deque<std::string> refIDs,
                              const std::vector<size_t> refLengths) const;
    size_t processReadRecords(const std::vector<SamRecord> &readRecords,
                              const std::deque<std::string> &referenceIDs, auto &splitsOut,
                              [[maybe_unused]] auto &multiSplitsOut) const;

    std::optional<SplitRecords> constructSplitRecords(const SamRecord &readRecord) const;
    std::optional<SplitRecords> constructSplitRecords(
        const std::vector<SamRecord> &readRecords) const;
    std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords> getSplitRecords(
        const std::vector<SamRecord> &readRecords,
        const std::deque<std::string> &referenceIDs) const;

    void mergeOutputFiles(const ChunkedOutTmpDirs &tmpDirs, const DetectOutput &output) const;
    void writeSamFile(auto &samOut, const std::vector<SamRecord> &splitRecords) const;

    void writeTranscriptCountsFile(const fs::path &transcriptCountsFilePath,
                                   const TranscriptCounts &transcriptCounts) const;
    void writeReadCountsSummaryFile(const Result &results, const std::string &sampleName,
                                    const fs::path &statsFilePath) const;

    ChunkedOutTmpDirs prepareTmpOutputDirs(const fs::path &tmpOutDir) const;
};

}  // namespace detect
}  // namespace pipelines
