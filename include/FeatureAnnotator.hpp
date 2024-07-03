#pragma once

// Boost
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard
#include <algorithm>
#include <optional>
#include <unordered_map>

// Classes
#include "DataTypes.hpp"
#include "FeatureParser.hpp"
#include "IITree.hpp"

namespace Annotation {

namespace uuids = boost::uuids;

using FeatureTreeMap = std::unordered_map<std::string, IITree<int, dtp::Feature>>;
class FeatureAnnotator {
   public:
    FeatureAnnotator(fs::path featureFilePath, const std::vector<std::string>& includedFeatures,
                     const std::string& featureIDFlag);
    FeatureAnnotator(fs::path featureFilePath, const std::vector<std::string>& includedFeatures);
    explicit FeatureAnnotator(const dtp::FeatureMap& featureMap);

    ~FeatureAnnotator() = default;

    class Results;

    size_t featureCount() const;
    std::string insert(const dtp::GenomicRegion& region);
    std::string mergeInsert(const dtp::GenomicRegion& region);

    std::vector<dtp::Feature> overlappingFeatures(const dtp::GenomicRegion& region);
    Results overlappingFeatureIterator(const dtp::GenomicRegion& region);
    std::optional<dtp::Feature> getBestOverlappingFeature(const dtp::GenomicRegion& region);

   private:
    FeatureTreeMap featureTreeMap;

    FeatureTreeMap buildFeatureTreeMap(const fs::path& featureFilePath,
                                       const std::vector<std::string>& includedFeatures,
                                       const std::optional<std::string>& featureIDFlag);
    FeatureTreeMap buildFeatureTreeMap(const dtp::FeatureMap& featureMap);
};

class FeatureAnnotator::Results {
   public:
    Results(const IITree<int, dtp::Feature>& tree, const std::vector<size_t>& indices,
            const std::optional<dtp::Strand>& strand);

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