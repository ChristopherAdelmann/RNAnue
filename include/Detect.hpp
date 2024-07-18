#pragma once

// openMP
#include <omp.h>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>

// Standard
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// Classes
#include "CooptimalPairwiseAligner.hpp"
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"
#include "EvaluatedSplitRecords.hpp"
#include "FeatureAnnotator.hpp"
#include "Logger.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

using namespace dtp;
class Detect {
   public:
    explicit Detect(po::variables_map params);
    ~Detect() = default;

    void start(pt::ptree sample);

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

    po::variables_map params;

    size_t minReadLength;
    double minComplementarity;
    double minFraction;
    int minMapQuality;
    bool excludeSoftClipping;
    bool filterSplicedReads;
    int spliceFilterTolerance;
    Annotation::Orientation annotationOrientation;
    Annotation::FeatureAnnotator featureAnnotator;

    Results iterateSortedMappingsFile(const std::string &mappingsInPath,
                                      const std::string &splitsPath,
                                      const std::string &multSplitsPath,
                                      const fs::path &unassignedSingletonRecordsOutPath);
    size_t processReadRecords(const std::vector<SamRecord> &readRecords,
                              const std::deque<std::string> &referenceIDs, auto &splitsOut,
                              [[maybe_unused]] auto &multiSplitsOut);

    std::optional<SplitRecords> constructSplitRecords(const SamRecord &readRecord);
    std::optional<SplitRecords> constructSplitRecords(const std::vector<SamRecord> &readRecords);
    std::optional<EvaluatedSplitRecords> getSplitRecords(
        const std::vector<SamRecord> &readRecords, const std::deque<std::string> &referenceIDs);

    void mergeResults(Results &transcriptCounts, const Results &newTranscriptCounts) const;
    void writeSamFile(auto &samOut, const std::vector<SamRecord> &splitRecords);

    void writeTranscriptCountsFile(const fs::path &transcriptCountsFilePath,
                                   const TranscriptCounts &transcriptCounts) const;
    void writeStatisticsFile(const Results &results, const std::string &sampleName,
                             const fs::path &statsFilePath) const;

    std::vector<ChunkedInOutFilePaths> prepareInputOutputFiles(const fs::path &mappingsFilePath,
                                                               const fs::path &splitsFilePath,
                                                               const int mappingRecordsCount);
    std::vector<fs::path> splitMappingsFile(const fs::path &mappingsFilePath,
                                            const fs::path &tmpInPath, const int entries);
};