#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>

// seqan3
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/views/chunk.hpp>

// Internal
#include "AnalyzeData.hpp"
#include "AnalyzeParameters.hpp"
#include "AnalyzeSample.hpp"
#include "FeatureAnnotator.hpp"
#include "InteractionCluster.hpp"

namespace math = boost::math;  // NOLINT

using seqan3::operator""_tag;

namespace pipelines::analyze {

class Analyze {
   public:
    explicit Analyze(AnalyzeParameters params)
        : parameters(params), featureAnnotator(params.featuresInPath, params.featureTypes) {};

    void process(const AnalyzeData &data);

   private:
    AnalyzeParameters parameters;
    annotation::FeatureAnnotator featureAnnotator;

    void processSample(AnalyzeSample sample);

    void assignTranscriptsToClusters(std::vector<InteractionCluster> &clusters,
                                     const std::deque<std::string> &referenceIDs,
                                     const AnalyzeSample &sample);
    static void assignAnnotatedContiguousFragmentCountsToTranscripts(
        const fs::path &contiguousTranscriptCountsInPath,
        std::unordered_map<std::string, size_t> &transcriptCounts);

    static auto getTranscriptProbabilities(
        const std::unordered_map<std::string, size_t> &transcriptCounts,
        size_t totalTranscriptCount) -> std::unordered_map<std::string, double>;

    static void assignPValuesToClusters(
        std::vector<InteractionCluster> &clusters,
        const std::unordered_map<std::string, size_t> &transcriptCounts,
        size_t totalTranscriptCount);
    static void assignPAdjustedValuesToClusters(std::vector<InteractionCluster> &clusters);
    void writeInteractionsToFile(const std::vector<InteractionCluster> &mergedClusters,
                                 const std::string &sampleName, const fs::path &clusterOutPath,
                                 const std::deque<std::string> &referenceIDs) const;

    [[nodiscard]] static inline auto getReferenceID(const int32_t referenceIDIndex,
                                                    const std::deque<std::string> &referenceIDs)
        -> std::string {
        return referenceIDs[referenceIDIndex];
    }

    static void writeBEDLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                                   const std::deque<std::string> &referenceIDs,
                                   const std::string &color, std::ofstream &bedOut);

    static void writeBEDArcLineToFile(const InteractionCluster &cluster,
                                      const std::string &clusterID,
                                      const std::deque<std::string> &referenceIDs,
                                      std::ofstream &bedOut);

    static void writeInteractionLineToFile(const InteractionCluster &cluster,
                                           const std::string &clusterID,
                                           const std::deque<std::string> &referenceIDs,
                                           std::ofstream &interactionOut);

    void assignNonAnnotatedContiguousToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath, annotation::FeatureAnnotator &featureAnnotator,
        std::unordered_map<std::string, size_t> &transcriptCounts);
    static auto parseSampleFragmentCount(const fs::path &sampleCountsInPath) -> size_t;
};

}  // namespace pipelines::analyze
