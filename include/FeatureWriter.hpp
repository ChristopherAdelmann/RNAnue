#pragma once

// Standard
#include <fstream>

// Classes
#include "DataTypes.hpp"
#include "FeatureAnnotator.hpp"
#include "FeatureParser.hpp"

namespace FeatureWriter {
void write(const Annotation::FeatureTreeMap &featureTreeMap, const std::string &outputPath,
           const Annotation::FileType::Value fileType);
};