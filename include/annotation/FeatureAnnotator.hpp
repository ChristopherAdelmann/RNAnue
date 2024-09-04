#pragma once

// Boost
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard
#include <deque>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Classes
#include "DataTypes.hpp"
#include "FeatureParser.hpp"
#include "IITree.hpp"
#include "Orientation.hpp"

namespace Annotation {
using FeatureTreeMap = std::unordered_map<std::string, IITree<int, dtp::Feature>>;
namespace fs = std::filesystem;

class FeatureAnnotator {
 public:
  FeatureAnnotator(fs::path featureFilePath,
                   const std::unordered_set<std::string>& includedFeatures,
                   const std::string& featureIDFlag);
  FeatureAnnotator(fs::path featureFilePath,
                   const std::unordered_set<std::string>& includedFeatures);

  explicit FeatureAnnotator(const dtp::FeatureMap& featureMap);
  explicit FeatureAnnotator() = default;

  ~FeatureAnnotator() = default;

  class Results;

  struct MergeInsertResult {
    std::string featureID;
    std::vector<std::string> mergedFeatureIDs;
  };

  size_t featureCount() const;
  std::string insert(const dtp::GenomicRegion& region);
  MergeInsertResult mergeInsert(const dtp::GenomicRegion& region, const int graceDistance);

  std::vector<dtp::Feature> overlappingFeatures(const dtp::GenomicRegion& region,
                                                const Orientation orientation);
  Results overlappingFeatureIterator(const dtp::GenomicRegion& region,
                                     const Orientation orientation) const;
  std::optional<dtp::Feature> getBestOverlappingFeature(const dtp::GenomicRegion& region,
                                                        const Orientation orientation) const;
  std::optional<dtp::Feature> getBestOverlappingFeature(const dtp::SamRecord& record,
                                                        const std::deque<std::string>& referenceIDs,
                                                        const Orientation orientation) const;
  const FeatureTreeMap& getFeatureTreeMap() const;

 private:
  FeatureTreeMap featureTreeMap;

  FeatureTreeMap buildFeatureTreeMap(const fs::path& featureFilePath,
                                     const std::vector<std::string>& includedFeatures,
                                     const std::optional<std::string>& featureIDFlag);
  FeatureTreeMap buildFeatureTreeMap(const fs::path& featureFilePath,
                                     const std::unordered_set<std::string>& includedFeatures,
                                     const std::optional<std::string>& featureIDFlag);
  FeatureTreeMap buildFeatureTreeMap(const dtp::FeatureMap& featureMap);
};

class FeatureAnnotator::Results {
 public:
  Results(const IITree<int, dtp::Feature>& tree, const std::vector<size_t>& indices,
          const std::optional<dtp::Strand>& strand);

  Results() = delete;

  struct Iterator;

  [[nodiscard]] Iterator begin() const;
  [[nodiscard]] Iterator end() const;

 private:
  const IITree<int, dtp::Feature>& tree;
  std::vector<size_t> indices;
  std::optional<dtp::Strand> strand;
};

struct FeatureAnnotator::Results::Iterator {
  using value_type = dtp::Feature;
  using difference_type = std::ptrdiff_t;
  using reference = const dtp::Feature&;
  using pointer = const dtp::Feature*;
  using iterator_category = std::forward_iterator_tag;

  explicit Iterator(const IITree<int, dtp::Feature>* tree, const std::vector<size_t>& indices,
                    size_t index, const std::optional<dtp::Strand>& strand);

  [[nodiscard]] reference operator*() const;
  [[nodiscard]] pointer operator->() const;

  Iterator& operator++();
  Iterator operator++(int);

  friend bool operator==(const Iterator& a, const Iterator& b) {
    return a.current_index == b.current_index;
  }

  friend bool operator!=(const Iterator& a, const Iterator& b) {
    return a.current_index != b.current_index;
  }

 private:
  const IITree<int, dtp::Feature>* tree;
  std::vector<size_t> indices;
  size_t current_index;
  std::optional<dtp::Strand> strand;
};

}  // namespace Annotation
