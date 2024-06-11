#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// htslib
#include <htslib/hts.h>

// segemehl
extern "C" {
#include <segemehl.h>
}

// Standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <optional>
#include <regex>
#include <sstream>
#include <utility>

// Class
#include "Logger.hpp"
#include "Utility.hpp"

// Samtools derived elements
extern "C" {
typedef enum {
    Coordinate,
    QueryName,
    TagCoordinate,
    TagQueryName,
    MinHash,
    TemplateCoordinate
} SamOrder;
int bam_sort_core_ext(SamOrder sam_order, char *sort_tag, int minimiser_kmer, bool try_rev,
                      bool no_squash, const char *fn, const char *prefix, const char *fnout,
                      const char *modeout, size_t _max_mem, int n_threads, const htsFormat *in_fmt,
                      const htsFormat *out_fmt, char *arg_list, int no_pg, int write_index);
}

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

class Align {
   public:
    explicit Align(po::variables_map params);
    ~Align() = default;

    void start(pt::ptree sample);

   private:
    po::variables_map params;
    fs::path indexPath;

    void buildIndex();
    void alignReads(const std::string &query, const std::string &mate, const std::string &matched);
    void sortAlignmentsByQueryName(const std::string &alignmentsPath,
                                   const std::string &sortedAlignmentsPath);
};
