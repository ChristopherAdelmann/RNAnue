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
#include "FeatureWriter.hpp"
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
    double maxComplementarityScore;
    double minHybridizationEnergy;

    static std::optional<Segment> fromSamRecord(const dtp::SamRecord &record);

    dtp::GenomicRegion toGenomicRegion(const std::deque<std::string> &referenceIDs) const;

    dtp::Feature toFeature(const std::deque<std::string> &referenceIDs,
                           const std::string &featureID,
                           const std::string &featureType = "transcript") const;

    void merge(const Segment &other);
};

struct InteractionCluster {
    std::pair<Segment, Segment> segments;
    std::vector<double> complementarityScores;
    std::vector<double> hybridizationEnergies;
    std::optional<std::pair<std::string, std::string>> transcriptIDs = std::nullopt;
    std::optional<double> pValue = std::nullopt;
    std::optional<double> pAdj = std::nullopt;
    int count{1};

    static std::optional<InteractionCluster> fromSegments(const Segment &segment1,
                                                          const Segment &segment2);

    bool operator<(const InteractionCluster &a) const;

    bool overlaps(const InteractionCluster &other, const int graceDistance) const;

    void merge(const InteractionCluster &other);

    double complementarityStatistics() const;

    double hybridizationEnergyStatistics() const;

   private:
    InteractionCluster(std::pair<Segment, Segment> segments,
                       const std::vector<double> &complementarityScores,
                       const std::vector<double> &hybridizationEnergies);

    static std::optional<std::pair<Segment, Segment>> getSortedElements(const Segment &segment1,
                                                                        const Segment &segment2);
};

class Analyze {
   public:
    explicit Analyze(po::variables_map params);
    ~Analyze() = default;
    void start(pt::ptree sample);

   private:
    po::variables_map params;
    Annotation::FeatureAnnotator featureAnnotator;

    void iterateSplitsFile(const fs::path &splitsInPath, const fs::path &unassignedSingletonsInPath,
                           const fs::path &fragmentCountsInPath, const fs::path &clusterOutPath,
                           const fs::path &supplementaryFeaturesOutPath,
                           const fs::path &clusterTranscriptCountsOutPath);
    void mergeOverlappingClusters(std::vector<InteractionCluster> &clusters);
    void assignClustersToTranscripts(std::vector<InteractionCluster> &clusters,
                                     const std::deque<std::string> &referenceIDs,
                                     const fs::path &unassignedSingletonsInPath,
                                     const fs::path &fragmentCountsInPath,
                                     const fs::path &supplementaryFeaturesOutPath,
                                     const fs::path &transcriptCountsOutPath);
    std::unordered_map<std::string, double> getTranscriptProbabilities(
        const std::unordered_map<std::string, size_t> &transcriptCounts,
        const size_t totalTranscriptCount);

    void assignPValuesToClusters(std::vector<InteractionCluster> &clusters,
                                 const std::unordered_map<std::string, size_t> &transcriptCounts,
                                 const size_t totalTranscriptCount);
    void assignPAdjustedValuesToClusters(std::vector<InteractionCluster> &clusters);
    void writeInteractionsToFile(const std::vector<InteractionCluster> &mergedClusters,
                                 const fs::path &clusterOutPath,
                                 const std::deque<std::string> &referenceIDs);

    inline std::string getReferenceID(const int32_t referenceIDIndex,
                                      const std::deque<std::string> &referenceIDs) const {
        return referenceIDs[referenceIDIndex];
    }

    void writeBEDLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                            const std::deque<std::string> &referenceIDs, const std::string &color,
                            std::ofstream &bedOut);

    void writeBEDArcLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                               const std::deque<std::string> &referenceIDs,
                               const std::string &color, std::ofstream &bedOut);

    void writeInteractionLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                                    const std::deque<std::string> &referenceIDs,
                                    std::ofstream &interactionOut);

    void assignNonAnnotatedSingletonsToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
        std::unordered_map<std::string, size_t> &transcriptCounts);
    size_t parseSampleFragmentCount(const fs::path &sampleCountsInPath);
};