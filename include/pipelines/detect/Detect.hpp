#pragma once
// openMP
#include <omp.h>

// Standard
#include <algorithm>
#include <deque>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>

// Boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// Classes
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"
#include "DetectData.hpp"
#include "DetectParameters.hpp"
#include "DetectSample.hpp"
#include "FeatureAnnotator.hpp"
#include "Logger.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsEvaluator.hpp"
#include "Utility.hpp"

using namespace dtp;

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
    struct ChunkedInOutFilePaths {
        fs::path mappingsInPath;
        fs::path splitsOutPath;
        fs::path multSplitsOutPath;
        fs::path unassignedSingletonRecordsOutPath;
    };
    using TranscriptCounts = std::unordered_map<std::string, size_t>;

    struct Results {
        TranscriptCounts transcriptCounts;
        size_t splitFragmentsCount{0};
        size_t singletonFragmentsCount{0};
    };

    DetectParameters params;

    Annotation::FeatureAnnotator featureAnnotator;
    SplitRecordsEvaluator splitRecordsEvaluator;

    const SplitRecordsEvaluationParameters::ParameterVariant getSplitRecordsEvaluatorParameters(
        const DetectParameters &params) const;

    const std::deque<std::string> &getReferenceIDs(const fs::path &mappingsInPath) const;

    void processSample(const DetectSample &sample) const;

    Results iterateSortedMappingsFile(const std::string &mappingsInPath,
                                      const std::string &splitsPath,
                                      const std::string &multiSplitsPath,
                                      const fs::path &unassignedSingletonRecordsOutPath) const;
    size_t processReadRecords(const std::vector<SamRecord> &readRecords,
                              const std::deque<std::string> &referenceIDs, auto &splitsOut,
                              [[maybe_unused]] auto &multiSplitsOut) const;

    std::optional<SplitRecords> constructSplitRecords(const SamRecord &readRecord) const;
    std::optional<SplitRecords> constructSplitRecords(
        const std::vector<SamRecord> &readRecords) const;
    std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords> getSplitRecords(
        const std::vector<SamRecord> &readRecords,
        const std::deque<std::string> &referenceIDs) const;

    void mergeResults(Results &transcriptCounts, const Results &newTranscriptCounts) const;
    void writeSamFile(auto &samOut, const std::vector<SamRecord> &splitRecords) const;

    void writeTranscriptCountsFile(const fs::path &transcriptCountsFilePath,
                                   const TranscriptCounts &transcriptCounts) const;
    void writeStatisticsFile(const Results &results, const std::string &sampleName,
                             const fs::path &statsFilePath) const;

    std::vector<ChunkedInOutFilePaths> prepareInputOutputFiles(const fs::path &mappingsFileInPath,
                                                               const fs::path &splitsFileOutPath,
                                                               const int mappingRecordsCount) const;
    std::vector<fs::path> splitMappingsFile(const fs::path &mappingsFilePath,
                                            const fs::path &tmpInPath, const int entries) const;
};

}  // namespace detect
}  // namespace pipelines
