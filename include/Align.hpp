// boost
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex>
#include <utility>
#include <filesystem>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

#include <bitset>

#include "Logger.hpp"

namespace pt = boost::property_tree;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace seqan3::literals;

typedef std::pair<uint32_t, uint32_t>
    ReadPos;
typedef std::pair<uint64_t,uint64_t> GenomePos;
typedef std::vector<seqan3::cigar> CigarSplt;

typedef std::vector<
    std::tuple<
        std::string,
        seqan3::sam_flag,
        std::optional<int32_t>,
        std::optional<int32_t>,
        std::vector<seqan3::cigar>,
        seqan3::dna5_vector,
        seqan3::sam_tag_dictionary>>
    Splts;

class Align {
    private:
        po::variables_map params;
        std::string segemehlSysCall;
        std::string index;

        void sortAlignments(std::string alignmentsPath);

    public:
        // constructor
        Align(po::variables_map params);
        Align();

        void buildIndex();
        void alignReads(std::string query, std::string matched);
        void detSplits(std::string matched, std::string splits);

        double complementarity(seqan3::dna5_vector rna1, seqan3::dna5_vector rna2);
        double hybridize(seqan3::dna5_vector rna1, seqan3::dna5_vector rna2);

        void processSplits(auto &splitrecords, auto &splitsfile);

        seqan3::dna5 string2dna5(std::string rna); 

        void constructIndex();
        void start(pt::ptree sample);
};
