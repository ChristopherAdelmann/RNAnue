#include "ParameterOptions.hpp"

#include <cstddef>

#include "boost/program_options/options_description.hpp"

namespace pi = constants::pipelines;

po::options_description ParameterOptions::getSubcallOptions() {
    po::options_description subcall("Subcall");
    subcall.add_options()("subcall", po::value<std::string>(), pi::SUBCALL_DESCRIPTION.c_str());

    return subcall;
}

po::options_description ParameterOptions::getGeneralOptions() {
    po::options_description general(pi::GENERAL_DESCRIPTION.c_str());
    general.add_options()("trtms,t", po::value<std::string>(),
                          "folder containing the raw reads of the treatments including replicates "
                          "located within sub-folders (condition) (required)");
    general.add_options()("ctrls,s", po::value<std::string>(),
                          "folder containing the the raw reads of the controls including "
                          "replicates located within sub-folders (condition)");
    general.add_options()("outdir,o", po::value<std::string>(),
                          "(output) folder in which the results are stored (required)");
    general.add_options()("loglevel", po::value<std::string>()->default_value("info"),
                          "log level [info, warning, error] (default: info)");
    general.add_options()("threads,p", po::value<int>()->default_value(1),
                          "max number of threads to be used (default: 1)");
    general.add_options()("features,f", po::value<std::string>(),
                          "annotation/features in GFF/GTF format (required)");
    general.add_options()(
        "featuretypes", po::value<std::string>()->default_value(std::string{"transcript"}),
        "feature types to be considered for the analysis, can be specified "
        "as --featuretypes 'gene,rRNA' comma seperated values (default: transcript)");
    general.add_options()(
        "orientation",
        po::value<Annotation::Orientation>()->default_value(Annotation::Orientation::BOTH),
        "orientation of the features to consider in relation to reads [same, "
        "opposite, both] (default: both)");
    general.add_options()("chunksize", po::value<int>()->default_value(100000),
                          "number of reads processed per chunk in parallel (default: 100000)");

    return general;
}

po::options_description ParameterOptions::getPreprocessOptions() {
    po::options_description preprocess("Preprocess Pipeline");
    preprocess.add_options()(pi::PREPROCESS.c_str(), po::bool_switch()->default_value(true),
                             "whether to include preprocessing of the raw reads in the "
                             "workflow of RNAnue (default: true)");
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
    preprocess.add_options()("minovltrim", po::value<size_t>()->default_value(5),
                             "minimum length of overlap between adapter and read (default: 5)");
    preprocess.add_options()(
        "minqual,q", po::value<size_t>()->default_value(20),
        "lower limit for the mean quality (Phred Quality Score) of the reads (default: 20)");
    preprocess.add_options()("minlen,l", po::value<size_t>()->default_value(15),
                             "minimum length of the reads (default: 15)");
    preprocess.add_options()(
        "wqual", po::value<size_t>()->default_value(20),
        "minimum mean quality for each window (Phred Quality Score) (default: 20)");
    preprocess.add_options()("wtrim", po::value<size_t>()->default_value(0),
                             "window size for quality trimming from 3' end. Selecting '0' will not "
                             "apply quality trimming (default: 0)");
    preprocess.add_options()("minovl", po::value<size_t>()->default_value(5),
                             "minimal overlap to merge paired-end reads (default: 5)");
    preprocess.add_options()(
        "mmerge", po::value<double>()->default_value(0.05),
        "rate of mismatches allowed when merging paired end reads (default: 0.05)");

    return preprocess;
}

po::options_description ParameterOptions::getAlignOptions() {
    po::options_description align("Align Pipeline");
    align.add_options()("dbref", po::value<std::string>(), "reference genome (.fasta) (required)");
    align.add_options()("accuracy", po::value<size_t>()->default_value(90),
                        "minimum percentage of read matches (default: 90, range: 0-100)");
    align.add_options()("minfragsco", po::value<size_t>()->default_value(18),
                        "minimum score of a spliced fragment (default: 18)");
    align.add_options()("minfraglen", po::value<size_t>()->default_value(20),
                        "minimum length of a spliced fragment (default: 20)");
    align.add_options()("minsplicecov", po::value<size_t>()->default_value(80),
                        "minimum coverage for spliced transcripts (default: 80, range:0-100)");

    return align;
}

po::options_description ParameterOptions::getDetectOptions() {
    po::options_description detect("Detect Pipeline");
    detect.add_options()("mapqmin", po::value<size_t>()->default_value(10),
                         "minimum quality of the alignments (default: 10)");
    detect.add_options()("cmplmin", po::value<double>()->default_value(0.0),
                         "complementarity cutoff for split reads (default: 0.0; range: 0.0-1.0)");
    detect.add_options()("sitelenratio", po::value<double>()->default_value(0.1),
                         "aligned portion of the read (default: 0.1, range: 0.0-1.0)");
    detect.add_options()(
        "nrgmax", po::value<double>()->default_value(0.0),
        "hybridization energy cutoff for split reads (default: 0.0, range: >=0.0)");
    detect.add_options()("exclclipping", po::bool_switch()->default_value(false),
                         "exclude soft clipping from the alignments (default: false)");
    detect.add_options()("splicing", po::bool_switch()->default_value(false),
                         "splicing events are removed in the detection of split reads");
    detect.add_options()("splicingtolerance", po::value<int>()->default_value(5),
                         "tolerance for splicing events (default: 5)");

    return detect;
}

po::options_description ParameterOptions::getAnalyzeOptions() {
    po::options_description analysis("Analyze Pipeline");
    analysis.add_options()("clustdist", po::value<int>()->default_value(0),
                           "threshold distance at which two clusters are merged into a single "
                           "combined cluster, default is to only merge overlapping or blunt "
                           "ended clusters (default: 0)");
    analysis.add_options()("padj", po::value<double>()->default_value(1.0),
                           "p-value threshold for outputting an interaction (default: 1.0)");
    analysis.add_options()("mincount", po::value<size_t>()->default_value(1),
                           "minimum number of reads assigned to an interaction (default: 1)");

    return analysis;
}

po::options_description ParameterOptions::getOtherOptions() {
    po::options_description other("Other");
    other.add_options()("version,v", "display the version number");
    other.add_options()("help,h", "display this help message");
    other.add_options()("config,c", po::value<std::string>(),
                        "configuration file that contains the parameters");

    return other;
}
