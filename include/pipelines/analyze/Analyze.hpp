#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <fstream>

// Standard
#include <algorithm>
#include <filesystem>
#include <optional>
#include <random>
#include <ranges>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

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
#include "AnalyzeData.hpp"
#include "AnalyzeParameters.hpp"
#include "AnalyzeSample.hpp"
#include "Constants.hpp"
#include "CustomSamTags.hpp"
#include "DataTypes.hpp"
#include "FeatureAnnotator.hpp"
#include "FeatureWriter.hpp"
#include "InteractionCluster.hpp"
#include "Logger.hpp"
#include "Segment.hpp"
#include "Utility.hpp"

namespace math = boost::math;  // NOLINT

using namespace dtp;
using seqan3::operator""_tag;

namespace pipelines {
namespace analyze {

class Analyze {
   public:
    explicit Analyze(AnalyzeParameters params)
        : parameters(params), featureAnnotator(params.featuresInPath, params.featureTypes) {};
    ~Analyze() = default;

    void process(const AnalyzeData &data);

   private:
    AnalyzeParameters parameters;
    Annotation::FeatureAnnotator featureAnnotator;

    void processSample(AnalyzeSample sample);

    void mergeOverlappingClusters(std::vector<InteractionCluster> &clusters);

    void assignClustersToTranscripts(std::vector<InteractionCluster> &clusters,
                                     const std::deque<std::string> &referenceIDs,
                                     const AnalyzeSample &sample);
    void assignAnnotatedContiguousFragmentCountsToTranscripts(
        const fs::path &contiguousTranscriptCountsInPath,
        std::unordered_map<std::string, size_t> &transcriptCounts);

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
                               const std::deque<std::string> &referenceIDs, std::ofstream &bedOut);

    void writeInteractionLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                                    const std::deque<std::string> &referenceIDs,
                                    std::ofstream &interactionOut);

    void assignNonAnnotatedContiguousToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
        std::unordered_map<std::string, size_t> &transcriptCounts);
    size_t parseSampleFragmentCount(const fs::path &sampleCountsInPath);
};

}  // namespace analyze
}  // namespace pipelines
