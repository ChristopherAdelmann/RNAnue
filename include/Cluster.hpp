#pragma once

// OpenMP
#include <omp.h>

// Standard
#include <algorithm>
#include <filesystem>
#include <random>
#include <ranges>

// Boost
#include <boost/filesystem.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// seqan3
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/views/chunk.hpp>

// Classes
#include "Constants.hpp"
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"
#include "FeatureAnnotator.hpp"
#include "Logger.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;
namespace math = boost::math;

using namespace dtp;
using seqan3::operator""_tag;

struct Segment {
    std::string recordID;
    int32_t referenceIDIndex;
    dtp::Strand strand;
    int32_t start;
    int32_t end;
    double complementarityScore;
    double hybridizationEnergy;

    static std::optional<Segment> fromSamRecord(const dtp::SamRecord &record);

    dtp::GenomicRegion toGenomicRegion(const std::deque<std::string> &referenceIDs) const;

    dtp::Feature toFeature(const std::deque<std::string> &referenceIDs,
                           const std::string &featureID,
                           const std::string &featureType = "transcript") const;

    void merge(const Segment &other);
};

struct ReadCluster {
    std::pair<Segment, Segment> segments;
    std::vector<double> complementarityScores;
    std::vector<double> hybridizationEnergies;
    std::optional<std::pair<std::string, std::string>> transcriptIDs = std::nullopt;
    std::optional<double> pValue = std::nullopt;
    std::optional<double> pAdj = std::nullopt;
    int count{1};

    static std::optional<ReadCluster> fromSegments(const Segment &segment1,
                                                   const Segment &segment2);

    bool operator<(const ReadCluster &a) const;

    bool overlaps(const ReadCluster &other, const int graceDistance) const;

    void merge(const ReadCluster &other);

    double complementarityStatistics() const;

    double hybridizationEnergyStatistics() const;

   private:
    explicit ReadCluster(std::pair<Segment, Segment> segments,
                         std::vector<double> complementarityScores,
                         std::vector<double> hybridizationEnergies);

    static std::optional<std::pair<Segment, Segment>> getSortedElements(const Segment &segment1,
                                                                        const Segment &segment2);
};

class Cluster {
   public:
    explicit Cluster(po::variables_map params);
    ~Cluster() = default;
    void start(pt::ptree sample);

   private:
    po::variables_map params;
    Annotation::FeatureAnnotator featureAnnotator;

    void iterateSplitsFile(const fs::path &splitsInPath, const fs::path &unassignedSingletonsInPath,
                           const fs::path &fragmentCountsInPath, const fs::path &clusterOutPath,
                           const fs::path &supplementaryFeaturesOutPath,
                           const fs::path &clusterTranscriptCountsOutPath);
    void mergeOverlappingClusters(std::vector<ReadCluster> &clusters);
    void assignClustersToTranscripts(std::vector<ReadCluster> &clusters,
                                     const std::deque<std::string> &referenceIDs,
                                     const fs::path &unassignedSingletonsInPath,
                                     const fs::path &fragmentCountsInPath,
                                     const fs::path &supplementaryFeaturesOutPath,
                                     const fs::path &transcriptCountsOutPath);
    std::unordered_map<std::string, double> getTranscriptProbabilities(
        const std::unordered_map<std::string, size_t> &transcriptCounts,
        const size_t totalTranscriptCount);

    void assignPValuesToClusters(std::vector<ReadCluster> &clusters,
                                 const std::unordered_map<std::string, size_t> &transcriptCounts,
                                 const size_t totalTranscriptCount);
    void assignPAdjustedValuesToClusters(std::vector<ReadCluster> &clusters);
    void writeClustersToFile(const std::vector<ReadCluster> &mergedClusters,
                             const fs::path &clusterOutPath,
                             const std::deque<std::string> &referenceIDs);
    void assignNonAnnotatedSingletonsToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
        std::unordered_map<std::string, size_t> &transcriptCounts);
    size_t parseSampleFragmentCount(const fs::path &sampleCountsInPath);
};
