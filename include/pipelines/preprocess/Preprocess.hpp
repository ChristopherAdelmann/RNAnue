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
#include "Constants.hpp"
#include "Logger.hpp"
#include "PairedRecordMerger.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"
#include "TrimConfig.hpp"
#include "Utility.hpp"
#include "VariantOverload.hpp"
#include "boost/filesystem/path.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/io/sequence_file/input.hpp"

using seqan3::operator""_dna5;

namespace pipelines {
namespace preprocess {
class Preprocess {
   public:
    explicit Preprocess(const PreprocessParameters &params);
    ~Preprocess() = default;

    void process(const PreprocessData &data) const;

   private:
    PreprocessParameters parameters;

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

    bool passesFilters(const auto &record) const;

    void processSample(const PreprocessSampleType &sample) const;

    void processSingleEnd(const PreprocessSampleSingle &sample) const;
    void processPairedEnd(const PreprocessSamplePaired &sample) const;

    void processSingleEndRecordChunk(SingleEndFastqChunk &chunk,
                                     const std::vector<Adapter> &adapters5,
                                     const std::vector<Adapter> &adapters3) const;
    void processSingleEndFileInChunks(std::string const &recInPath, std::string recOutPath,
                                      const std::vector<Adapter> &adapters5,
                                      const std::vector<Adapter> &adapters3, size_t chunkSize,
                                      size_t numThreads) const;

    void processPairedEndRecordChunk(PairedEndFastqChunk &chunk,
                                     const std::vector<Adapter> &adapters5f,
                                     const std::vector<Adapter> &adapters3f,
                                     const std::vector<Adapter> &adapters5r,
                                     const std::vector<Adapter> &adapters3r) const;
    void processPairedEndFileInChunks(
        std::string const &recFwdInPath, std::string const &recRevInPath,
        std::string const &mergedOutPath, std::string const &snglFwdOutPath,
        std::string const &snglRevOutPath, const std::vector<Adapter> &adapters5f,
        const std::vector<Adapter> &adapters3f, const std::vector<Adapter> &adapters5r,
        const std::vector<Adapter> &adapters3r, size_t chunkSize, size_t numThreads) const;
};

}  // namespace preprocess
}  // namespace pipelines
