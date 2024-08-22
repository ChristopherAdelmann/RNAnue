#pragma once
// Standard
#include <string>

namespace constants {
namespace pipelines {
const std::string PREPROCESS = "preprocess";
const std::string ALIGN = "align";
const std::string DETECT = "detect";
const std::string ANALYZE = "analyze";
const std::string COMPLETE = "complete";

const std::string GENERAL_DESCRIPTION =
    "RNAnue efficient data analysis for RNA–RNA interactomics.\nRun RNAnue with the subcall "
    "\"complete\" to execute all pipeline steps.\n\nMinimum call: RNAnue complete -t "
    "<treatment-dir> "
    "-o "
    "<output-dir> -f <feature-gff-file> --dbref <reference-genome-file>\nOr run RNAnue with a "
    "config "
    "file: RNAnue complete -c <config-file>\n\nGeneral Options";
const std::string SUBCALL_DESCRIPTION =
    "The subcall to execute. The following subcalls are available: preprocess, align, detect, "
    "analyze, complete.";

const std::string PROCESSING_TREATMENT_MESSAGE = "Processing treatment data";
const std::string PROCESSING_CONTROL_MESSAGE = "Processing control data";
}  // namespace pipelines
}  // namespace constants
