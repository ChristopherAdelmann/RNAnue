#ifndef RNANUE_ALIGN_HPP
#define RNANUE_ALIGN_HPP

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

// Class
#include "Utility.hpp"

// Standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bitset>
#include <filesystem>
#include <iostream>
#include <regex>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <sstream>
#include <utility>

#include "Logger.hpp"

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace seqan3::literals;

typedef std::pair<uint32_t, uint32_t> ReadPos;
typedef std::pair<uint64_t, uint64_t> GenomePos;
typedef std::vector<seqan3::cigar> CigarSplt;

typedef std::vector<
    std::tuple<std::string, seqan3::sam_flag, std::optional<int32_t>, std::optional<int32_t>,
               std::vector<seqan3::cigar>, seqan3::dna5_vector, seqan3::sam_tag_dictionary>>
    Splits;
#include <bitset>
#include <filesystem>

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

class Align {
   public:
    // constructor/destructor
    Align(po::variables_map params);
    ~Align() = default;

    // alignment
    void buildIndex();
    void alignReads(std::string query, std::string mate, std::string matched);
    void start(pt::ptree sample);
    seqan3::dna5 string2dna5(std::string rna);

   private:
    po::variables_map params;
    std::string segemehlSysCall;
    std::string index;

    void sortAlignments(std::string alignmentsPath);
};

#endif  // RNANUE_ALIGN_HPP
