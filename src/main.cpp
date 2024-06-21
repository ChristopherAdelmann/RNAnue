
// Boost
#include <boost/program_options.hpp>

// Standard
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <array>
#include <iostream>
#include <string>

// Classes
#include "Base.hpp"
#include "Closing.hpp"
#include "Config.hpp"
#include "Constants.hpp"

namespace po = boost::program_options;

void showVersion(std::ostream& _str) {
    _str << "RNAnue v" << RNAnue_VERSION_MAJOR;
    _str << "." << RNAnue_VERSION_MINOR << ".";
    _str << RNAnue_VERSION_PATCH << " - ";
    _str << "Detect RNA-RNA interactions ";
    _str << "from Direct-Duplex-Detection (DDD) data.";
    _str << std::endl;
}

// A helper function to simplify the main part.
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

namespace std {
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vec) {
    for (auto item : vec) {
        os << item << " ";
    }
    return os;
}
}  // namespace std

void handler(int sig) {
    std::array<void*, 10> array;
    size_t size = backtrace(array.data(), 10);

    // print out all the frames to stderr
    std::cerr << "Error: signal " << sig << ":" << std::endl;
    backtrace_symbols_fd(array.data(), size, STDERR_FILENO);
    exit(1);
}

int main(int argc, char* argv[]) {
    signal(SIGSEGV, handler);  // register handler to catch segmentation faults

    try {
        using namespace constants::pipelines;

        std::string readType;
        std::string configFile;

        po::options_description general("General");
        general.add_options()("readtype,r", po::value<std::string>(&readType)->default_value("SE"),
                              "single-end (=SE) or paired-end (=PE) reads");
        general.add_options()(
            "trtms,t", po::value<std::string>(),
            "folder containing the raw reads of the treatments including replicates "
            "located within subfolders (condition)");
        general.add_options()("ctrls,s", po::value<std::string>(),
                              "folder containing the the raw reads of the controls including "
                              "replicates located within subfolders (condition)");
        general.add_options()("outdir,o", po::value<std::string>(),
                              "(output) folder in which the results are stored");
        general.add_options()("loglevel", po::value<std::string>()->default_value("info"),
                              "log level (info, warning, error)");
        general.add_options()("threads,p", po::value<int>()->default_value(1),
                              "max number of threads");
        general.add_options()("features", po::value<std::string>(),
                              "annotation/features in .GFF3 format");
        general.add_options()("featuretypes",
                              po::value<std::vector<std::string>>()->multitoken()->default_value(
                                  std::vector<std::string>{"gene"}),
                              "feature types to be considered for the analysis, can be specified "
                              "as --featuretypes gene rRNA ... (default: gene)");
        general.add_options()(
            "mapqmin", po::value<int>()->default_value(20),
            "lower limit for the quality (Phred Quality Score) of the alignments");
        // general.add_options()("splicing", po::bool_switch->default_value(false),
        //                       "splicing events are considered in the detection of split
        //                       reads");

        po::options_description preprocess("Preprocessing");
        preprocess.add_options()(
            pi::PREPROCESS.c_str(), po::bool_switch()->default_value(true),
            "include preprocessing of the raw reads in the workflow of RNAnue");
        preprocess.add_options()("chunksize", po::value<int>()->default_value(100000),
                                 "number of reads to be processed in parallel (default: 100000)");
        preprocess.add_options()(
            "trimpolyg", po::bool_switch()->default_value(false),
            "trim high quality polyG tails from the reads. Applicable for Illumina "
            "NextSeq reads. (default: false)");
        preprocess.add_options()(
            "adpt5f", po::value<std::string>()->default_value(""),
            "single sequence or file [.fasta] of the adapter sequences to be removed "
            "from the 5' end (forward read if PE)");
        preprocess.add_options()(
            "adpt5r", po::value<std::string>()->default_value(""),
            "single sequence or file [.fasta] of the adapter sequences to be removed "
            "from the 5' end of the reverse read (PE only)");
        preprocess.add_options()(
            "adpt3f", po::value<std::string>()->default_value(""),
            "single sequence or file [.fasta] of the adapter sequences to be removed "
            "from the 3' end (forward read if PE)");
        preprocess.add_options()(
            "adpt3r", po::value<std::string>()->default_value(""),
            "single sequence or file [.fasta] of the adapter sequences to be removed "
            "from the 3' end of the reverse read (PE only)");
        preprocess.add_options()(
            "mtrim", po::value<double>()->default_value(0.05),
            "rate of mismatches allowed when aligning adapters to sequences (default: 0.05)");
        preprocess.add_options()("mintrim", po::value<int>()->default_value(5),
                                 "minimum length of overlap between adapter and read (default: 5)");
        preprocess.add_options()(
            "minqual,q", po::value<int>()->default_value(20),
            "lower limit for the mean quality (Phred Quality Score) of the reads (default: 20)");
        preprocess.add_options()("minlen,l", po::value<std::size_t>()->default_value(15),
                                 "minimum length of the reads");
        preprocess.add_options()(
            "wqual", po::value<std::size_t>()->default_value(20),
            "minimum mean quality for each window (Phred Quality Score) (default: 20)");
        preprocess.add_options()(
            "wtrim", po::value<std::size_t>()->default_value(0),
            "window size for quality trimming from 3' end. Selecting '0' will not "
            "apply quality trimming (default: 0)");
        preprocess.add_options()("minovl", po::value<int>()->default_value(5),
                                 "minimal overlap to merge paired-end reads (default: 5)");
        preprocess.add_options()(
            "mmerge", po::value<double>()->default_value(0.05),
            "rate of mismatches allowed when merging paired end reads (default: 0.05)");

        po::options_description align("Alignment");
        align.add_options()("dbref", po::value<std::string>(), "reference genome (.fasta)")(
            "accuracy", po::value<int>()->default_value(90), "minimum percentage of read matches")(
            "minfragsco", po::value<int>()->default_value(18),
            "minimum score of a spliced fragment")("minfraglen",
                                                   po::value<int>()->default_value(20),
                                                   "minimum length of a spliced fragment")(
            "minsplicecov", po::value<int>()->default_value(80),
            "minimum coverage for spliced transcripts")(
            "exclclipping", po::bool_switch()->default_value(false),
            "exclude soft clipping from the alignments");  //(
        // "unprd", po::bool_switch()->default_value(false),
        // "only for paired-end reads: include unpaired reads")(
        // "unmrg", po::bool_switch()->default_value(false),
        // "only for paired-end reads: include unmerged reads");

        po::options_description detect("Split Read Calling");
        detect.add_options()("cmplmin", po::value<double>()->default_value(0.0),
                             "complementarity cutoff for split reads")(
            "sitelenratio", po::value<double>()->default_value(0.1),
            "aligned portion of the read (complementarity)")(
            "nrgmax", po::value<double>()->default_value(0),
            "hybridization energy cutoff for split reads");

        po::options_description clustering("Clustering");
        clustering.add_options()
            // ("clust", po::bool_switch()->default_value(true),
            //                          "include clustering of the split reads in the workflow of
            //                          RNAnue")
            ("clustdist", po::value<int>()->default_value(0), "minimum distance between clusters");

        po::options_description analysis("Analysis");

        po::options_description output("Output");
        output.add_options()(
            "stats", po::bool_switch()->default_value(false),
            "whether (=1) or not (=0) to (additionally) create statistics of the libraries")(
            "outcnt", po::bool_switch()->default_value(false),
            "whether (=1) or not (=0) to (additionally) save results as count table for DEA");
        //(
        // "outjgf", po::bool_switch()->default_value(false),
        // "whether (=1) or not (=0) to (additionally) save results as JSON Graph FILE (JGF)");

        po::options_description other("Other");
        other.add_options()("version,v", "display the version number")("help,h",
                                                                       "display this help message")(
            "config,c", po::value<std::string>(&configFile)->default_value("params.cfg"),
            "configuration file that contains the parameters");

        po::options_description subcall("Subcall");
        subcall.add_options()("subcall", po::value<std::string>(), SUBCALL_DESCRIPTION.c_str());

        po::options_description cmdlineOptions;
        cmdlineOptions.add(general)
            .add(preprocess)
            .add(align)
            .add(detect)
            .add(clustering)
            .add(analysis)
            .add(other)
            .add(subcall);

        po::options_description configFileOptions;
        configFileOptions.add(general)
            .add(preprocess)
            .add(align)
            .add(detect)
            .add(clustering)
            .add(analysis)
            .add(output);

        // translate all positional options into subcall options
        po::positional_options_description p;
        p.add("subcall", 1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).options(cmdlineOptions).positional(p).run(), vm);
        notify(vm);

        Logger::setLogLevel(vm["loglevel"].as<std::string>());

        // include parameters from the configfile if available
        std::ifstream ifs(configFile.c_str());

        Closing cl{};  // class that handles the closing remarks
        if (vm.count("help")) {
            std::cout << cmdlineOptions << std::endl;
            cl.printQuote(std::cout);
            return 0;
        }

        if (vm.count("version")) {
            showVersion(std::cout);
            cl.printQuote(std::cout);
            return 0;
        }

        if (!ifs) {
            Logger::log(LogLevel::ERROR, "Configuration file could not be opened!");
            return 1;
        } else {
            po::store(po::parse_config_file(ifs, configFileOptions), vm);
            notify(vm);
        }

        if (!vm.count("subcall")) {
            Logger::log(LogLevel::ERROR, "Please provide a subcall.");
            return 1;
        }

        // start execution
        showVersion(std::cout);
        Base bs(vm);  // controls all downstream processing

        cl.printQuote(std::cout);
    } catch (po::error& e) {
        std::cout << "please provide a correct function call" << std::endl;
        std::cerr << e.what() << " at line " << __LINE__ << std::endl;
        return 0;
    }
}
