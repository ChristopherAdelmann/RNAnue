#pragma once

// Standard
#include <algorithm>
#include <iostream>
#include <mutex>
#include <numeric>
#include <optional>
#include <ranges>
#include <thread>
#include <vector>

// boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan3
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/sequence_file/all.hpp>

// Class
#include "Logger.hpp"
#include "Utility.hpp"

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

using seqan3::operator""_dna5;
class Preprocess {
   public:
    explicit Preprocess(const po::variables_map &params);
    ~Preprocess() = default;

    void start(pt::ptree sample);

   private:
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

    struct TrimConfig {
       public:
        enum Mode { FIVE_PRIME, THREE_PRIME };

        static seqan3::align_cfg::method_global alignmentConfigFor(Mode mode);
    };

    struct Adapter {
        const seqan3::dna5_vector sequence;
        const double maxMissMatchFraction;
        const TrimConfig::Mode trimmingMode;
    };

    std::vector<Adapter> loadAdapters(const std::string &filenameOrSequence,
                                      const TrimConfig::Mode trimmingMode);

    bool passesFilters(const auto &record);

    template <typename record_type>
    void trimWindowedQuality(record_type &record);

    template <typename record_type>
    void trimAdapter(const Adapter &adapter, record_type &record);

    template <typename record_type>
    void trim3PolyG(record_type &record);

    template <typename record_type>
    std::optional<record_type> mergeRecordPair(const record_type &record1,
                                               const record_type &record2);

    template <typename record_type, typename result_type>
    record_type constructMergedRecord(const record_type &record1, const record_type &record2,
                                      const seqan3::alignment_result<result_type> &alignmentResult);

    void processSingleEnd(pt::ptree sample);
    void processPairedEnd(pt::ptree sample);

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