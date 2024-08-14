#pragma once

// Standard
#include <algorithm>
#include <iostream>
#include <mutex>
#include <numeric>
#include <optional>
#include <ranges>
#include <thread>
#include <variant>
#include <vector>

// boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan3
#include <seqan3/io/sequence_file/all.hpp>

// Class
#include "Adapter.hpp"
#include "Logger.hpp"
#include "PairedRecordMerger.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"
#include "TrimConfig.hpp"
#include "Utility.hpp"

namespace pt = boost::property_tree;

using seqan3::operator""_dna5;

template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

namespace pipelines {
namespace preprocess {
class Preprocess {
   public:
    explicit Preprocess(const PreprocessParameters &params);
    ~Preprocess() = default;

    void start(pt::ptree sample);
    void start(const PreprocessData &data);

   private:
    PreprocessParameters parameters;
    std::string readType;

    bool trimPolyG;

    std::string adpt5f;
    std::string adpt5r;
    std::string adpt3f;
    std::string adpt3r;
    double missMatchRateTrim;
    int minOverlapTrim;

    int minMeanPhread;
    std::size_t minLen;

    std::size_t minWindowPhread;
    std::size_t windowTrimSize;

    int minOverlapMerge;
    double missMatchRateMerge;

    size_t threadCount;
    size_t chunkSize;

    struct SingleEndFastqChunk {
        std::vector<seqan3::sequence_file_input<>::record_type> records;
    };

    struct PairedEndFastqChunk {
        std::vector<seqan3::sequence_file_input<>::record_type> recordsFwd;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsRev;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsMergedRes;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsSnglFwdRes;
        std::vector<seqan3::sequence_file_input<>::record_type> recordsSnglRevRes;
    };

    bool passesFilters(const auto &record);

    void processSample(const PreprocessSampleType &sample);

    void processSingleEnd(const PreprocessSampleSingle &sample);
    void processPairedEnd(const PreprocessSamplePaired &sample);

    void processSingleEndRecordChunk(SingleEndFastqChunk &chunk,
                                     const std::vector<Adapter> &adapters5,
                                     const std::vector<Adapter> &adapters3);
    void processSingleEndFileInChunks(std::string const &recInPath, std::string recOutPath,
                                      const std::vector<Adapter> &adapters5,
                                      const std::vector<Adapter> &adapters3, size_t chunkSize,
                                      size_t numThreads);

    void processPairedEndRecordChunk(PairedEndFastqChunk &chunk,
                                     const std::vector<Adapter> &adapters5f,
                                     const std::vector<Adapter> &adapters3f,
                                     const std::vector<Adapter> &adapters5r,
                                     const std::vector<Adapter> &adapters3r);
    void processPairedEndFileInChunks(
        std::string const &recFwdInPath, std::string const &recRevInPath,
        std::string const &mergedOutPath, std::string const &snglFwdOutPath,
        std::string const &snglRevOutPath, const std::vector<Adapter> &adapters5f,
        const std::vector<Adapter> &adapters3f, const std::vector<Adapter> &adapters5r,
        const std::vector<Adapter> &adapters3r, size_t chunkSize, size_t numThreads);
};
}  // namespace preprocess
}  // namespace pipelines
