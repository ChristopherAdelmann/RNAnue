#ifndef SEQRICKSHAW_H
#define SEQRICKSHAW_H

#include "Helper.hpp"
#include "Logger.hpp"

// openMP
#include <omp.h>

#include <bitset>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <ranges>
#include <thread>
#include <mutex>
#include <optional>

// boost
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>

// seqan
#include <seqan3/io/sequence_file/all.hpp>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/all.hpp>

#include <seqan3/io/exception.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/slice.hpp>

#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <range/v3/all.hpp>
#include <range/v3/view/transform.hpp>

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

typedef std::tuple<std::string,int,std::pair<int,int>,std::size_t> State;
typedef std::vector<State> States;

typedef std::map<std::pair<int,char>,std::tuple<int,int,int,int>> LookupTable;

typedef std::map<std::pair<std::string,std::string>, LookupTable> Adapters;

typedef std::vector<fs::path> PathVector;
typedef std::pair<PathVector,PathVector> PathVectorPair;

//typedef std::tuple<seqan3::field::id,seqan3::field::seq,seqan3::field:qual> 

//typedef std::map<std::pair<int,string::string>
//

using seqan3::operator""_dna5;
using seqan3::operator""_dna4;

class TrimConfig
{
public:
    enum Mode
    {
        FIVE_PRIME,
        THREE_PRIME
    };

    static seqan3::align_cfg::method_global alignmentConfigFor(Mode mode);
};
class SeqRickshaw {
    private:
        int modus; // defines the modus of the trimming procedure (0 == 5' trimming, 1 == 3' trimming, 2 == both)

        po::variables_map params;

        std::string adpt5prime;
        std::string adpt3prime;

        std::optional<Adapters> adpt5Table;
        std::optional<Adapters> adpt3Table;

        int phred; // the minimum
        int minLen;
        int wsize;

        std::string readtype;

        void setupLookupTables();

        std::chrono::duration<double> cumDuration;
        std::size_t sumReads;

        struct FastqChunk
        {
            std::vector<seqan3::sequence_file_input<>::record_type> records;
            // Add any other data you need from the FASTQ file.
        };

        void processChunk(FastqChunk const &chunk);
        void processFastqFileInChunks(std::string const &filename, size_t chunkSize, size_t numThreads);

        // Define a structure to hold the data of a single FASTQ chunk.

        /* the lookup table - smart transition table
         * (state,c) -> shift, state, readPos, match
         * */
        //        std::map<std::pair<int,char>,std::tuple<int,int,int,int>> lookup;

    public:
        SeqRickshaw(const po::variables_map &params);
        SeqRickshaw();

        std::map<std::pair<std::string,std::string>,LookupTable> calcLookupTable(std::string _type, std::string _path);
        LookupTable calcShift(auto &_records);

        std::vector<char> determineAlphabet(auto _sequence);
        std::size_t calcReadPos(auto& sequence, std::size_t& left, std::pair<std::size_t,std::size_t>& right);


        void smallestShift(std::string pattern, std::string suffix, int left);
        int transition(std::string pattern, std::string suffix, int readPos, std::size_t& left, std::pair<std::size_t,std::size_t>& right);
        std::pair<std::string, std::string> merging(auto fwd, auto rev, auto fwdQual, auto revQual);

        std::string longestCommonSubstr(std::string forward, std::string reverse);

        // #### 
        // split the readsfile specified in 'path' and returns a list of files
        std::vector<fs::path> splitReadsFile(fs::path path, std::string subfolder);
        std::vector<fs::path> splitOutputFile(std::vector<fs::path> inputfiles, std::string folder);
        void rmTmpOutDirs(fs::path outPath, std::vector<std::string> subfolders); // remove temporary output directories

        void distribute(std::pair<fs::path,fs::path> reads, fs::path merged, std::pair<fs::path,fs::path> unmerged, fs::path unpaired);
        void distribute(fs::path input, fs::path output);

        // helper 
        int addState(States &states, State state, States::size_type &size);
        int nextReadPos(std::string state, int currReadPos);
        // finds all occurrences of substring in string
        void findAllOcc(std::vector<std::size_t>& fnd, std::string str, std::string substr);
        void writeLookupTables(std::ofstream &os, Adapters &adptLookTbls);

        void preprocPattern();

		// perform window trimming
		std::size_t nibble(auto &seq, auto &qual, std::pair<std::size_t,std::size_t> &bnds);

		int trimmingWindow(auto seq);
        
        std::size_t boyermoore(auto& read, LookupTable tab, int patlen);
        std::pair<std::size_t,std::size_t> trimming(auto& read);
        void start(pt::ptree sample);

        bool passesFilters(auto &record);
        void trimAdapter(const seqan3::dna5_vector &adapterSequence, auto &recordSequence, TrimConfig::Mode trimmingMode);

        template <typename record_type>
        void trim3PolyG(record_type &record);

        template <typename record_type>
        std::optional<record_type> mergeRecordPair(record_type &record1, record_type &record2, int minOverlap);

        template <typename record_type>
        record_type constructMergedRecord(const record_type &record1, const record_type &record2, const int overlap);
};
#endif
