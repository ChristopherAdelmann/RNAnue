#pragma once

// Standard
#include <string>

namespace constants {
namespace pipelines {
constexpr std::string PREPROCESS = "preprocess";
constexpr std::string ALIGN = "align";
constexpr std::string DETECT = "detect";
constexpr std::string CLUSTER = "cluster";
constexpr std::string ANALYZE = "analyze";
constexpr std::string COMPLETE = "complete";

const std::string SUBCALL_DESCRIPTION =
    "The subcall to execute. The following subcalls are available: preprocess, align, detect, "
    "cluster, analyze, complete.The complete subcall will execute all subcalls in the correct "
    "order.";
}  // namespace pipelines
}  // namespace constants