#include "FeatureAnnotator.hpp"

using namespace Annotation;

FeatureAnnotator::FeatureAnnotator(fs::path featureFilePath,
                                   const std::vector<std::string> &includedFeatures,
                                   const std::string &featureIDFlag)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, featureIDFlag)) {}

FeatureAnnotator::FeatureAnnotator(fs::path featureFilePath,
                                   const std::vector<std::string> &includedFeatures)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, std::nullopt)) {}

FeatureAnnotator::FeatureAnnotator(const dtp::FeatureMap &featureMap)
    : featureTreeMap(buildFeatureTreeMap(featureMap)) {}

Annotation::FeatureTreeMap FeatureAnnotator::buildFeatureTreeMap(
    const dtp::FeatureMap &featureMap) {
    Annotation::FeatureTreeMap newFeatureTreeMap;
    newFeatureTreeMap.reserve(featureMap.size());

    for (const auto &[seqid, features] : featureMap) {
        IITree<int, dtp::Feature> tree;
        for (const auto &feature : features) {
            tree.add(feature.startPosition, feature.endPosition, feature);
        }
        tree.index();
        newFeatureTreeMap.emplace(seqid, std::move(tree));
    }

    return newFeatureTreeMap;
}

Annotation::FeatureTreeMap FeatureAnnotator::buildFeatureTreeMap(
    const fs::path &featureFilePath, const std::vector<std::string> &includedFeatures,
    const std::optional<std::string> &featureIDFlag) {
    Annotation::FeatureTreeMap newFeatureTreeMap;

    const std::unordered_set<std::string> includedFeaturesSet(includedFeatures.begin(),
                                                              includedFeatures.end());

    dtp::FeatureMap featureMap =
        FeatureParser(includedFeaturesSet, featureIDFlag).parse(featureFilePath);

    newFeatureTreeMap.reserve(featureMap.size());

    for (const auto &[seqid, features] : featureMap) {
        IITree<int, dtp::Feature> tree;
        for (const auto &feature : features) {
            tree.add(feature.startPosition, feature.endPosition, feature);
        }
        tree.index();
        newFeatureTreeMap.emplace(seqid, std::move(tree));
    }

    return newFeatureTreeMap;
}

std::vector<dtp::Feature> FeatureAnnotator::overlappingFeatures(const dtp::GenomicRegion &region) {
    std::vector<dtp::Feature> features;
    if (auto it = featureTreeMap.find(region.referenceID); it != featureTreeMap.end()) {
        std::vector<size_t> indices;
        it->second.overlap(region.startPosition, region.endPosition, indices);
        for (const auto &index : indices) {
            features.push_back(it->second.data(index));
        }
    }
    return features;
}

FeatureAnnotator::Results FeatureAnnotator::overlappingFeatureIterator(
    const dtp::GenomicRegion &region) {
    std::vector<size_t> indices;
    auto it = featureTreeMap.find(region.referenceID);
    if (it != featureTreeMap.end()) {
        it->second.overlap(region.startPosition, region.endPosition, indices);
    }
    return Results(it->second, std::move(indices));
}

std::optional<dtp::Feature> FeatureAnnotator::getBestOverlappingFeature(
    const dtp::GenomicRegion &region) {
    auto overlapSizeWithRegion = [region](const dtp::Feature &feature) -> size_t {
        return std::min(region.endPosition, feature.endPosition) -
               std::max(region.startPosition, feature.startPosition);
    };

    auto featureIterator = overlappingFeatureIterator(region);

    auto maxOverlapSizeElement = std::max_element(
        featureIterator.begin(), featureIterator.end(),
        [&overlapSizeWithRegion](const dtp::Feature &lhs, const dtp::Feature &rhs) {
            return std::invoke(overlapSizeWithRegion, lhs) <
                   std::invoke(overlapSizeWithRegion, rhs);
        });

    if (maxOverlapSizeElement != featureIterator.end()) {
        return *maxOverlapSizeElement;
    }

    return std::nullopt;
}

// Results and Iterator implementation
FeatureAnnotator::Results::Results(const IITree<int, dtp::Feature> &tree,
                                   const std::vector<size_t> &indices)
    : tree(tree), indices(indices) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::begin() const {
    return Iterator(&tree, indices, 0);
}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::end() const {
    return Iterator(&tree, indices, indices.size());
}

FeatureAnnotator::Results::Iterator::Iterator(const IITree<int, dtp::Feature> *tree,
                                              const std::vector<size_t> &indices, size_t index)
    : tree(tree), indices(indices), current_index(index) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator::reference
FeatureAnnotator::Results::Iterator::operator*() const {
    return tree->data(indices[current_index]);
}

[[nodiscard]] FeatureAnnotator::Results::Iterator::pointer
FeatureAnnotator::Results::Iterator::operator->() const {
    return &tree->data(indices[current_index]);
}

FeatureAnnotator::Results::Iterator &FeatureAnnotator::Results::Iterator::operator++() {
    ++current_index;
    return *this;
}

FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::Iterator::operator++(int) {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}