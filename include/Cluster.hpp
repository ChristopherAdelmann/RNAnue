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
#include "DataTypes.hpp"
#include "Logger.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

using namespace dtp;
struct Segment {
    std::string recordID;
    int32_t referenceIDIndex;
    uint32_t flag;
    uint32_t start;
    uint32_t end;

    static std::optional<Segment> fromSamRecord(const dtp::SamRecord &record) {
        if (!record.reference_position().has_value() || !record.reference_id().has_value()) {
            return std::nullopt;
        }

        const auto isReverseStrand =
            static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
        const uint32_t flag{isReverseStrand ? (u_int32_t)16 : (u_int32_t)0};

        const uint32_t start = record.reference_position().value();
        const uint32_t end = start + record.sequence().size() - 1;

        return Segment{record.id(), record.reference_id().value(), flag, start, end};
    }
};

struct ReadCluster {
    std::pair<Segment, Segment> elements;
    int count{1};

    ~ReadCluster() = default;

    bool operator<(const ReadCluster &a) const {
        return std::tie(elements.first.start, elements.second.start) <
               std::tie(a.elements.first.start, a.elements.second.start);
    }

    static std::optional<ReadCluster> fromSegments(const Segment &segment1,
                                                   const Segment &segment2) {
        const auto sortedElements = getSortedElements(segment1, segment2);
        if (!sortedElements.has_value()) {
            return std::nullopt;
        }
        return ReadCluster{sortedElements.value()};
    }

    bool overlaps(const ReadCluster &other, const int graceDistance) const {
        const bool firstOverlaps =
            elements.first.end + graceDistance >= other.elements.first.start &&
            elements.first.start <= other.elements.first.end + graceDistance;

        const bool secondOverlaps =
            elements.second.end + graceDistance >= other.elements.second.start &&
            elements.second.start <= other.elements.second.end + graceDistance;

        return elements.first.referenceIDIndex == other.elements.first.referenceIDIndex &&
               elements.second.referenceIDIndex == other.elements.second.referenceIDIndex &&
               firstOverlaps && secondOverlaps &&
               elements.first.flag == other.elements.first.flag &&
               elements.second.flag == other.elements.second.flag;
    }

    void merge(const ReadCluster &other) {
        elements.first.start = std::min(elements.first.start, other.elements.first.start);
        elements.first.end = std::max(elements.first.end, other.elements.first.end);
        elements.second.start = std::min(elements.second.start, other.elements.second.start);
        elements.second.end = std::max(elements.second.end, other.elements.second.end);
        count += other.count;
    }

   private:
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

    void iterate(std::string splits);
    void mergeOverlappingClusters(std::vector<ReadCluster> &clusters);
    void writeClustersToFile(const std::deque<std::string> &referenceIDs);
};
