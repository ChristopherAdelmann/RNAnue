#pragma once

// OpenMP
#include <omp.h>

// Standard
#include <algorithm>
#include <filesystem>
#include <ranges>

// Boost
#include <boost/filesystem.hpp>
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

using namespace dtp;
using seqan3::operator""_tag;

struct Segment {
    std::string recordID;
    int32_t referenceIDIndex;
    dtp::Strand strand;
    uint32_t start;
    uint32_t end;
    double complementarityScore;
    double hybridizationEnergy;

    static std::optional<Segment> fromSamRecord(const dtp::SamRecord &record) {
        if (!record.reference_position().has_value() || !record.reference_id().has_value()) {
            return std::nullopt;
        }

        const auto isReverseStrand =
            static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
        const dtp::Strand strand{isReverseStrand ? dtp::Strand::REVERSE : dtp::Strand::FORWARD};

        const uint32_t start = record.reference_position().value();
        const uint32_t end = start + record.sequence().size() - 1;

        const double hybridizationEnergy = record.tags().get<"XE"_tag>();
        const double complementarityScore = record.tags().get<"XC"_tag>();

        return Segment{record.id(),
                       record.reference_id().value(),
                       strand,
                       start,
                       end,
                       complementarityScore,
                       hybridizationEnergy};
    }

    dtp::GenomicRegion toGenomicRegion(const std::deque<std::string> &referenceIDs) const {
        return dtp::GenomicRegion{referenceIDs[referenceIDIndex], start, end};
    }

    dtp::Feature toFeature(const std::deque<std::string> &referenceIDs,
                           const std::string &featureID,
                           const std::string &featureType = "transcript") const {
        return dtp::Feature{
            referenceIDs[referenceIDIndex], featureType, start, end, strand, featureID};
    }

    void merge(const Segment &other) {
        start = std::min(start, other.start);
        end = std::max(end, other.end);
        complementarityScore = std::max(complementarityScore, other.complementarityScore);
        hybridizationEnergy = std::min(hybridizationEnergy, other.hybridizationEnergy);
    }
};

struct ReadCluster {
    std::pair<Segment, Segment> segments;
    std::vector<double> complementarityScores;
    std::vector<double> hybridizationEnergies;
    std::optional<std::pair<std::string, std::string>> transcriptIDs = std::nullopt;
    int count{1};

    bool operator<(const ReadCluster &a) const {
        return std::tie(segments.first.start, segments.second.start) <
               std::tie(a.segments.first.start, a.segments.second.start);
    }

    static std::optional<ReadCluster> fromSegments(const Segment &segment1,
                                                   const Segment &segment2) {
        const auto sortedElements = getSortedElements(segment1, segment2);
        if (!sortedElements.has_value()) {
            return std::nullopt;
        }
        assert(sortedElements.value().first.complementarityScore ==
               sortedElements.value().second.complementarityScore);
        assert(sortedElements.value().first.hybridizationEnergy ==
               sortedElements.value().second.hybridizationEnergy);
        return ReadCluster{sortedElements.value(),
                           {sortedElements->first.complementarityScore},
                           {sortedElements->first.hybridizationEnergy}};
    }

    bool overlaps(const ReadCluster &other, const int graceDistance) const {
        const bool firstOverlaps =
            segments.first.end + graceDistance >= other.segments.first.start &&
            segments.first.start <= other.segments.first.end + graceDistance;

        const bool secondOverlaps =
            segments.second.end + graceDistance >= other.segments.second.start &&
            segments.second.start <= other.segments.second.end + graceDistance;

        return segments.first.referenceIDIndex == other.segments.first.referenceIDIndex &&
               segments.second.referenceIDIndex == other.segments.second.referenceIDIndex &&
               firstOverlaps && secondOverlaps &&
               segments.first.strand == other.segments.first.strand &&
               segments.second.strand == other.segments.second.strand;
    }

    void merge(const ReadCluster &other) {
        segments.first.merge(other.segments.first);
        segments.second.merge(other.segments.second);
        complementarityScores.insert(complementarityScores.end(),
                                     other.complementarityScores.begin(),
                                     other.complementarityScores.end());
        hybridizationEnergies.insert(hybridizationEnergies.end(),
                                     other.hybridizationEnergies.begin(),
                                     other.hybridizationEnergies.end());
        count += other.count;
    }

    double complementarityStatistics() const {
        assert(!complementarityScores.empty());
        assert(segments.first.complementarityScore == segments.second.complementarityScore);
        return helper::calculateMedian(complementarityScores) / segments.first.complementarityScore;
    }

    double hybridizationEnergyStatistics() const {
        assert(!hybridizationEnergies.empty());
        assert(segments.first.hybridizationEnergy == segments.second.hybridizationEnergy);
        return std::sqrt(helper::calculateMedian(hybridizationEnergies) /
                         segments.first.hybridizationEnergy);
    }

   private:
    explicit ReadCluster(std::pair<Segment, Segment> segments,
                         std::vector<double> complementarityScores,
                         std::vector<double> hybridizationEnergies)
        : segments(segments),
          complementarityScores(complementarityScores),
          hybridizationEnergies(hybridizationEnergies) {}
    static std::optional<std::pair<Segment, Segment>> getSortedElements(const Segment &segment1,
                                                                        const Segment &segment2) {
        if (segment1.recordID != segment2.recordID) {
            Logger::log(LogLevel::WARNING, "Record IDs do not match: ", segment1.recordID, " vs. ",
                        segment2.recordID, ". Make sure reads are sorted by name.");
            return std::nullopt;
        }
        return segment1.start < segment2.start ? std::make_pair(segment1, segment2)
                                               : std::make_pair(segment2, segment1);
    }
};

class Cluster {
   public:
    explicit Cluster(po::variables_map params);
    ~Cluster() = default;
    void start(pt::ptree sample);

   private:
    po::variables_map params;
    std::vector<ReadCluster> result;
    Annotation::FeatureAnnotator featureAnnotator;

    void iterate(const std::string &splitsInPath, const fs::path &unassignedSingletonsInPath,
                 const fs::path &clusterOutPath, const fs::path &supplementaryFeaturesOutPath,
                 const fs::path &clusterTranscriptCountsOutPath);
    void mergeOverlappingClusters(std::vector<ReadCluster> &clusters);
    void assignClustersToTranscripts(const std::vector<ReadCluster> &clusters,
                                     const std::deque<std::string> &referenceIDs,
                                     const fs::path &unassignedSingletonsInPath,
                                     const fs::path &supplementaryFeaturesOutPath,
                                     const fs::path &transcriptCountsOutPath);
    void writeClustersToFile(const fs::path &clusterOutPath,
                             const std::deque<std::string> &referenceIDs);
    void assignUnassignedSingletonsToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath, Annotation::FeatureAnnotator &featureAnnotator,
        std::unordered_map<std::string, size_t> &transcriptCounts);
};
