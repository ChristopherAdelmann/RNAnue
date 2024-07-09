#pragma once

// Standard
#include <fstream>

// Classes
#include "DataTypes.hpp"
#include "FeatureAnnotator.hpp"
#include "FeatureParser.hpp"

namespace Annotation {
class FeatureWriter {
   public:
    FeatureWriter() = delete;
    ~FeatureWriter() = delete;
    static void write(const FeatureTreeMap &featureTreeMap, const std::string &outputPath,
                      const FileType::Value fileType);
};
}  // namespace Annotation