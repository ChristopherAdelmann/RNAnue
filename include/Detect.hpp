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
#include <format>
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
    };
    po::variables_map params;

    size_t minReadLength;
    bool excludeSoftClipping;
    double minComplementarity;
    double minFraction;

    void iterateSortedMappingsFile(const std::string &mappingsInPath, const std::string &splitsPath,
                                   const std::string &multSplitsPath);
    void processReadRecords(const std::vector<SamRecord> &readRecords, auto &splitsOut,
                            [[maybe_unused]] auto &multiSplitsOut);

    std::optional<SplitRecords> constructSplitRecords(const SamRecord &readRecord);
    std::optional<SplitRecords> constructSplitRecords(const std::vector<SamRecord> &readRecords);
    std::optional<EvaluatedSplitRecords> getSplitRecords(const std::vector<SamRecord> &readRecords);

    void writeSamFile(auto &samOut, const std::vector<SamRecord> &splitRecords);
    void createStatisticsFile(const fs::path &splitsFilePath, const fs::path &multSplitsFilePath,
                              const fs::path &statsFilePath) const;

    std::vector<ChunkedInOutFilePaths> prepareInputOutputFiles(const fs::path &mappingsFilePath,
                                                               const fs::path &splitsFilePath,
                                                               const int mappingRecordsCount);
    std::vector<fs::path> splitMappingsFile(const fs::path &mappingsFilePath,
                                            const fs::path &tmpInPath, const int entries) const;
};