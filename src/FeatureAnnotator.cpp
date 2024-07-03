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

    for (const auto &[referenceID, features] : featureMap) {
        IITree<int, dtp::Feature> tree;
        for (const auto &feature : features) {
            tree.add(feature.startPosition, feature.endPosition, feature);
        }
        tree.index();
        newFeatureTreeMap.emplace(referenceID, std::move(tree));
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

size_t FeatureAnnotator::featureCount() const {
    size_t count = 0;
    for (const auto &[_, tree] : featureTreeMap) {
        count += tree.size();
    }
    return count;
}

std::string FeatureAnnotator::insert(const dtp::GenomicRegion &region) {
    assert(region.strand.has_value() && "Strand must be specified for insertion");

    auto &tree = featureTreeMap[region.referenceID];
    const std::string uuid = uuids::to_string(uuids::random_generator()());
    tree.add(region.startPosition, region.endPosition,
             dtp::Feature{.referenceID = region.referenceID,
                          .type = "supplementary_feature",
                          .startPosition = region.startPosition,
                          .endPosition = region.endPosition,
                          .strand = *region.strand,
                          .id = uuid});
    tree.index();

    return uuid;
}

std::string FeatureAnnotator::mergeInsert(const dtp::GenomicRegion &region) {
    assert(region.strand.has_value() && "Strand must be specified for insertion");

    auto &tree = featureTreeMap[region.referenceID];

    std::vector<size_t> indices;
    tree.overlap(region.startPosition, region.endPosition, indices);

    std::erase_if(indices, [&region, &tree](size_t index) {
        return tree.data(index).strand != *region.strand;
    });

    if (indices.empty()) {
        return insert(region);
    }

    auto minStart = std::ranges::min_element(indices, [&tree](size_t lhs, size_t rhs) {
        return tree.data(lhs).startPosition < tree.data(rhs).startPosition;
    });

    auto maxEnd = std::ranges::max_element(indices, [&tree](size_t lhs, size_t rhs) {
        return tree.data(lhs).endPosition < tree.data(rhs).endPosition;
    });

    tree.data(*minStart).startPosition =
        std::min(tree.data(*minStart).startPosition, region.startPosition);
    tree.data(*minStart).endPosition = std::max(tree.data(*maxEnd).endPosition, region.endPosition);

    for (size_t i = indices.size() - 1; i > 0; --i) {
        if (i != *minStart) {
            tree.remove(indices[i]);
        }
    }

    tree.index();
    return tree.data(*minStart).id;
}

std::vector<dtp::Feature> FeatureAnnotator::overlappingFeatures(const dtp::GenomicRegion &region) {
    std::vector<dtp::Feature> features;
    if (auto it = featureTreeMap.find(region.referenceID); it != featureTreeMap.end()) {
        std::vector<size_t> indices;
        it->second.overlap(region.startPosition, region.endPosition, indices);
        for (const auto &index : indices) {
            if (region.strand == std::nullopt || it->second.data(index).strand == *region.strand) {
                features.push_back(it->second.data(index));
            }
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
    return Results(it->second, std::move(indices), region.strand);
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
                                   const std::vector<size_t> &indices,
                                   const std::optional<dtp::Strand> &strand)
    : tree(tree), indices(indices), strand(strand) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::begin() const {
    auto matchStrand = [&](size_t index) {
        return !strand.has_value() || tree.data(indices[index]).strand == strand;
    };
    size_t startIndex = 0;
    if (strand) {
        auto it = std::find_if(indices.begin(), indices.end(), matchStrand);
        startIndex = it != indices.end() ? std::distance(indices.begin(), it) : indices.size();
    }
    return Iterator(&tree, indices, startIndex, strand);
}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::end() const {
    return Iterator(&tree, indices, indices.size(), strand);
}

FeatureAnnotator::Results::Iterator::Iterator(const IITree<int, dtp::Feature> *tree,
                                              const std::vector<size_t> &indices, size_t index,
                                              const std::optional<dtp::Strand> &strand)
    : tree(tree), indices(indices), current_index(index), strand(strand) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator::reference
FeatureAnnotator::Results::Iterator::operator*() const {
    if (current_index >= indices.size()) throw std::out_of_range("Iterator out of range");
    return tree->data(indices[current_index]);
}

[[nodiscard]] FeatureAnnotator::Results::Iterator::pointer
FeatureAnnotator::Results::Iterator::operator->() const {
    if (current_index >= indices.size()) throw std::out_of_range("Iterator out of range");
    return &tree->data(indices[current_index]);
}

FeatureAnnotator::Results::Iterator &FeatureAnnotator::Results::Iterator::operator++() {
    auto matchStrand = [&](size_t index) {
        return !strand.has_value() || tree->data(indices[index]).strand == strand;
    };
    do {
        ++current_index;
    } while (current_index < indices.size() && !matchStrand(current_index));
    return *this;
}

FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::Iterator::operator++(int) {
    Iterator tmp = *this;
    ++(*this);  // Use pre-increment to simplify
    return tmp;
}