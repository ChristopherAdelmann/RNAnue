#pragma once

// Standard
#include <string>

namespace constants {
namespace pipelines {
const std::string PREPROCESS = "preprocess";
const std::string ALIGN = "align";
const std::string DETECT = "detect";
const std::string CLUSTER = "cluster";
const std::string ANALYZE = "analyze";
const std::string COMPLETE = "complete";

const std::string SUBCALL_DESCRIPTION =
    "The subcall to execute. The following subcalls are available: preprocess, align, detect, "
    "cluster, analyze, complete.The complete subcall will execute all subcalls in the correct "
    "order.";
}  // namespace pipelines
}  // namespace constants