#include "FeatureAnnotator.hpp"

#include <ios>
#include <iostream>
#include <unordered_set>

// Feature Annotator
using namespace Annotation;
FeatureAnnotator::FeatureAnnotator(fs::path featureFilePath,
                                   const std::vector<std::string> &includedFeatures,
                                   const std::string &featureIDFlag)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, featureIDFlag)) {}

FeatureAnnotator::FeatureAnnotator(fs::path featureFilePath,
                                   const std::vector<std::string> &includedFeatures)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, std::nullopt)) {}

FeatureAnnotator::FeatureAnnotator(fs::path featureFilePath,
                                   const std::unordered_set<std::string> &includedFeatures)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, std::nullopt)) {}

FeatureAnnotator::FeatureAnnotator(const dtp::FeatureMap &featureMap)
    : featureTreeMap(buildFeatureTreeMap(featureMap)) {}

FeatureTreeMap FeatureAnnotator::buildFeatureTreeMap(const dtp::FeatureMap &featureMap) {
    FeatureTreeMap newFeatureTreeMap;
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

FeatureTreeMap FeatureAnnotator::buildFeatureTreeMap(
    const fs::path &featureFilePath, const std::vector<std::string> &includedFeatures,
    const std::optional<std::string> &featureIDFlag) {
    FeatureTreeMap newFeatureTreeMap;

    std::unordered_set<std::string> uniqueIncludedFeatures;

    for (const std::string &feature : includedFeatures) {
        std::stringstream ss(feature);
        std::string str;
        while (getline(ss, str, ',')) {
            uniqueIncludedFeatures.insert(str);
        }
    }

    dtp::FeatureMap featureMap =
        FeatureParser(uniqueIncludedFeatures, featureIDFlag).parse(featureFilePath);

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

FeatureTreeMap FeatureAnnotator::buildFeatureTreeMap(
    const fs::path &featureFilePath, const std::unordered_set<std::string> &includedFeatures,
    const std::optional<std::string> &featureIDFlag) {
    FeatureTreeMap newFeatureTreeMap;

    dtp::FeatureMap featureMap =
        FeatureParser(includedFeatures, featureIDFlag).parse(featureFilePath);

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

const FeatureTreeMap &FeatureAnnotator::getFeatureTreeMap() const { return featureTreeMap; }

std::string FeatureAnnotator::insert(const dtp::GenomicRegion &region) {
    assert(region.strand.has_value() && "Strand must be specified for insertion");
    namespace uuids = boost::uuids;

    auto &tree = featureTreeMap[region.referenceID];
    const std::string uuid = uuids::to_string(uuids::random_generator()());
    tree.add(region.startPosition, region.endPosition,
             dtp::Feature{.referenceID = region.referenceID,
                          .type = "supplementary_feature",
                          .startPosition = region.startPosition,
                          .endPosition = region.endPosition,
                          .strand = *region.strand,
                          .id = uuid,
                          .groupID = std::nullopt});
    tree.index();

    return uuid;
}

FeatureAnnotator::MergeInsertResult FeatureAnnotator::mergeInsert(const dtp::GenomicRegion &region,
                                                                  const int graceDistance) {
    assert(region.strand.has_value() && "Strand must be specified for insertion");

    auto &tree = featureTreeMap[region.referenceID];

    std::vector<size_t> indices;
    // Overlap with grace distance and blunt ends (+/- 1)
    tree.overlap(region.startPosition - graceDistance - 1, region.endPosition + graceDistance + 1,
                 indices);

    std::erase_if(indices, [&region, &tree](size_t index) {
        return tree.data(index).strand != *region.strand;
    });

    if (indices.empty()) {
        return {insert(region), {}};
    }

    auto minStartIndex = *std::ranges::min_element(indices, [&tree](size_t lhs, size_t rhs) {
        return tree.data(lhs).startPosition < tree.data(rhs).startPosition;
    });

    auto maxEndIndex = *std::ranges::max_element(indices, [&tree](size_t lhs, size_t rhs) {
        return tree.data(lhs).endPosition < tree.data(rhs).endPosition;
    });

    dtp::Feature &minStartFeature = tree.data(minStartIndex);
    minStartFeature.startPosition = std::min(minStartFeature.startPosition, region.startPosition);
    minStartFeature.endPosition = std::max(tree.data(maxEndIndex).endPosition, region.endPosition);
    tree.setStart(minStartIndex, minStartFeature.startPosition);
    tree.setEnd(minStartIndex, minStartFeature.endPosition);

    std::vector<std::string> mergedFeatureIDs;

    for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
        if (*it != minStartIndex) {
            mergedFeatureIDs.push_back(tree.data(*it).id);
            tree.remove(*it);
        }
    }

    tree.index();
    return {.featureID = minStartFeature.id, .mergedFeatureIDs = mergedFeatureIDs};
}

std::vector<dtp::Feature> FeatureAnnotator::overlappingFeatures(const dtp::GenomicRegion &region,
                                                                const Orientation orientation) {
    std::vector<dtp::Feature> features;
    auto it = featureTreeMap.find(region.referenceID);
    if (it != featureTreeMap.end()) {
        std::vector<size_t> indices;
        it->second.overlap(region.startPosition, region.endPosition, indices);

        features.reserve(indices.size());

        for (const auto &index : indices) {
            const auto &feature = it->second.data(index);
            bool shouldAddFeature = false;

            if (region.strand == std::nullopt) {
                shouldAddFeature = true;
            } else if ((orientation == Orientation::BOTH) ||
                       (orientation == Orientation::OPPOSITE && feature.strand != *region.strand) ||
                       (orientation == Orientation::SAME && feature.strand == *region.strand)) {
                shouldAddFeature = true;
            }

            if (shouldAddFeature) {
                features.push_back(feature);
            }
        }
    }

    return features;
}

FeatureAnnotator::Results FeatureAnnotator::overlappingFeatureIterator(
    const dtp::GenomicRegion &region, const Orientation orientation) const {
    std::vector<size_t> indices;

    auto it = featureTreeMap.find(region.referenceID);

    if (it == featureTreeMap.end()) [[unlikely]] {
        return Results(it->second, std::move(indices), std::nullopt);
    }

    it->second.overlap(region.startPosition, region.endPosition, indices);

    std::optional<dtp::Strand> strand = std::nullopt;
    if (orientation == Orientation::SAME) {
        strand = region.strand;
    } else if (orientation == Orientation::OPPOSITE && region.strand.has_value()) {
        strand = !*region.strand;
    }

    return Results(it->second, std::move(indices), strand);
}

std::optional<dtp::Feature> FeatureAnnotator::getBestOverlappingFeature(
    const dtp::GenomicRegion &region, const Orientation orientation) const {
    auto overlapSizeWithRegion = [region](const dtp::Feature &feature) -> size_t {
        return std::min(region.endPosition, feature.endPosition) -
               std::max(region.startPosition, feature.startPosition);
    };

    auto featureIterator = overlappingFeatureIterator(region, orientation);

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

std::optional<dtp::Feature> FeatureAnnotator::getBestOverlappingFeature(
    const dtp::SamRecord &record, const std::deque<std::string> &referenceIDs,
    const Orientation orientation) const {
    auto region = dtp::GenomicRegion::fromSamRecord(record, referenceIDs);
    if (!region.has_value()) {
        return std::nullopt;
    }

    return getBestOverlappingFeature(region.value(), orientation);
}

// Results and Iterator implementation
FeatureAnnotator::Results::Results(const IITree<int, dtp::Feature> &tree,
                                   const std::vector<size_t> &indices,
                                   const std::optional<dtp::Strand> &strand)
    : tree(tree), indices(indices), strand(strand) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::begin() const {
    size_t startIndex = 0;
    if (strand.has_value()) {
        // Find the first index with the specified strand
        auto it = std::find_if(indices.begin(), indices.end(),
                               [&](size_t index) { return tree.data(index).strand == *strand; });
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
        if (strand.has_value()) {
            return tree->data(indices[index]).strand == *strand;
        }
        return true;
    };
    do {
        ++current_index;
    } while (current_index < indices.size() && !matchStrand(current_index));
    return *this;
}

FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::Iterator::operator++(int) {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}
