#ifndef RNANUE_CLUSTERING_HPP
#define RNANUE_CLUSTERING_HPP

#include <omp.h>

#include <algorithm>
#include <filesystem>

// boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/views/chunk.hpp>

#include "Logger.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

struct Segment {
    std::string refid;
    uint32_t flag;
    uint32_t start;
    uint32_t end;

    Segment() : refid(""), flag(0), start(0), end(0){};
    Segment(std::string _refid, uint32_t _flag, unsigned int _start, unsigned int _end)
        : refid(_refid), flag(_flag), start(_start), end(_end){};
};

struct ReadCluster {
    std::vector<Segment> elements;
    int count;

    ReadCluster() : elements({}), count(1){};
    bool operator<(const ReadCluster &a) const { return elements[0].start < a.elements[0].start; }
};

class Cluster {
   private:
    po::variables_map params;
    std::vector<ReadCluster> result;

   public:
    Cluster(po::variables_map params);

    void iterate(std::string splits);
    void overlaps(std::vector<ReadCluster> &clusterlist);
    bool startPosCmp(ReadCluster &a, ReadCluster &b);
    void start(pt::ptree sample);
    void sumup();
};

#endif  // RNANUE_CLUSTERING_HPP
