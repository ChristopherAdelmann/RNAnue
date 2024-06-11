#pragma once

// Boost
#include <boost/filesystem.hpp>

// Standard
#include <cstddef>
#include <unordered_map>

// Classes
#include "DataTypes.hpp"
#include "FeatureParser.hpp"
#include "IITree.hpp"

namespace Annotation {

using FeatureTreeMap = std::unordered_map<std::string, IITree<int, dtp::Feature>>;

class FeatureAnnotator {
   public:
    explicit FeatureAnnotator(fs::path featureFilePath,
                              const std::unordered_set<std::string>& includedFeatures,
                              const std::optional<std::string>& featureIDFlag);
    ~FeatureAnnotator() = default;

    class Results;

    std::vector<dtp::Feature> overlappingFeatures(const dtp::GenomicRegion& region);
    std::optional<Results> overlappingFeatureIt(const dtp::GenomicRegion& region);

   private:
    FeatureTreeMap featureTreeMap;

    FeatureTreeMap buildFeatureTreeMap(const fs::path& featureFilePath,
                                       const std::unordered_set<std::string>& includedFeatures,
                                       const std::optional<std::string>& featureIDFlag);
};

class FeatureAnnotator::Results {
   public:
    Results(const IITree<int, dtp::Feature>& tree, const std::vector<size_t>& indices);

    struct Iterator;

    [[nodiscard]] Iterator begin() const;
    [[nodiscard]] Iterator end() const;

    [[nodiscard]] Iterator cbegin() const;
    [[nodiscard]] Iterator cend() const;

   private:
    const IITree<int, dtp::Feature>& tree;
    std::vector<size_t> indices;
};

struct FeatureAnnotator::Results::Iterator {
    using value_type = dtp::Feature;
    using difference_type = std::ptrdiff_t;
    using reference = const dtp::Feature&;
    using pointer = const dtp::Feature*;
    using iterator_category = std::forward_iterator_tag;

    Iterator(const IITree<int, dtp::Feature>& tree, const std::vector<size_t>& indices,
             size_t index);

    [[nodiscard]] reference operator*() const;
    [[nodiscard]] pointer operator->() const;

    Iterator& operator++();
    Iterator operator++(int);

    friend bool operator==(const Iterator& a, const Iterator& b) {
        return a.current_index_ == b.current_index_;
    }

    friend bool operator!=(const Iterator& a, const Iterator& b) {
        return a.current_index_ != b.current_index_;
    }

   private:
    const IITree<int, dtp::Feature>& tree_;
    const std::vector<size_t> indices_;
    size_t current_index_;
};

}  // namespace Annotation