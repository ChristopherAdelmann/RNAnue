#pragma once

// Standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <optional>
#include <vector>

// Boost
#include <boost/program_options.hpp>

// htslib
#include <htslib/hts.h>

// segemehl
extern "C" {
#include <segemehl.h>
}  // segemehl

// Class
#include "AlignData.hpp"
#include "AlignParameters.hpp"
#include "AlignSample.hpp"

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

namespace pipelines {
namespace align {

class Align {
   public:
    explicit Align(AlignParameters params) : parameters(params) {};
    ~Align() = default;

    void process(const AlignData &data);

   private:
    AlignParameters parameters;
    fs::path indexPath;

    void processSample(const AlignSampleType &sample);

    void processSingleEnd(const AlignSampleSingle &sample);
    void processMergedPairedEnd(const AlignSampleMergedPaired &sample);

    std::optional<fs::path> findIndex(const fs::path &referenceGenomePath) const;

    void buildIndex();

    std::vector<std::string> getGeneralAlignmentArgs() const;
    void alignReads(const std::string &query, const std::string &mate,
                    const std::string &matched) const;
    void alignSingleReads(const fs::path &queryFastqInPath,
                          const fs::path &alignmentsFastqOutPath) const;
    void alignPairedReads(const fs::path &queryForwardFastqInPath,
                          const fs::path &queryReverseFastqInPath,
                          const fs::path &alignmentsFastqOutPath) const;

    void sortAlignmentsByQueryName(const fs::path &alignmentsPath,
                                   const fs::path &sortedAlignmentsPath) const;

    static auto convertToCStrings(std::vector<std::string> &args) -> std::vector<char *>;
};

}  // namespace align
}  // namespace pipelines
