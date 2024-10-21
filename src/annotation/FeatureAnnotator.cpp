#include "FeatureAnnotator.hpp"

// Standard
#include <ranges>
#include <unordered_set>

// Internal
#include "FeatureParser.hpp"
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"

namespace annotation {

FeatureAnnotator::FeatureAnnotator(fs::path &featureFilePath,
                                   const std::unordered_set<std::string> &includedFeatures,
                                   const std::string &featureIDFlag)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, featureIDFlag)) {}

FeatureAnnotator::FeatureAnnotator(fs::path &featureFilePath,
                                   const std::unordered_set<std::string> &includedFeatures)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, includedFeatures, std::nullopt)) {}

FeatureAnnotator::FeatureAnnotator(const dataTypes::FeatureMap &featureMap)
    : featureTreeMap(buildFeatureTreeMap(featureMap)) {}

auto FeatureAnnotator::buildFeatureTreeMap(const dataTypes::FeatureMap &featureMap)
    -> FeatureTreeMap {
    FeatureTreeMap newFeatureTreeMap;
    newFeatureTreeMap.reserve(featureMap.size());

    for (const auto &[referenceID, features] : featureMap) {
        IITree<int, dataTypes::Feature> tree;
        for (const auto &feature : features) {
            tree.add(feature.startPosition, feature.endPosition, feature);
        }
        tree.index();
        newFeatureTreeMap.emplace(referenceID, std::move(tree));
    }

    return newFeatureTreeMap;
}

auto FeatureAnnotator::buildFeatureTreeMap(
    const fs::path &featureFilePath, const std::unordered_set<std::string> &includedFeatures,
    const std::optional<std::string> &featureIDFlag) -> FeatureTreeMap {
    FeatureTreeMap newFeatureTreeMap;

    dataTypes::FeatureMap featureMap =
        FeatureParser(includedFeatures, featureIDFlag).parse(featureFilePath);

    newFeatureTreeMap.reserve(featureMap.size());

    for (const auto &[seqid, features] : featureMap) {
        IITree<int, dataTypes::Feature> tree;
        for (const auto &feature : features) {
            tree.add(feature.startPosition, feature.endPosition, feature);
        }
        tree.index();
        newFeatureTreeMap.emplace(seqid, std::move(tree));
    }

    return newFeatureTreeMap;
}

auto FeatureAnnotator::featureCount() const -> size_t {
    size_t count = 0;
    for (const auto &[_, tree] : featureTreeMap) {
        count += tree.size();
    }
    return count;
}

auto FeatureAnnotator::insert(const dataTypes::GenomicRegion &region) -> std::string {
    assert(region.strand.has_value() && "Strand must be specified for insertion");
    namespace uuids = boost::uuids;

    auto &tree = featureTreeMap[region.referenceID];
    const std::string uuid = uuids::to_string(uuids::random_generator()());
    tree.add(region.startPosition, region.endPosition,
             dataTypes::Feature{.referenceID = region.referenceID,
                                .type = "supplementary_feature",
                                .startPosition = region.startPosition,
                                .endPosition = region.endPosition,
                                .strand = *region.strand,
                                .id = uuid,
                                .groupID = std::nullopt,
                                .geneName = std::nullopt});
    tree.index();

    return uuid;
}

auto FeatureAnnotator::mergeInsert(const dataTypes::GenomicRegion &region,
                                   const int graceDistance) -> FeatureAnnotator::MergeInsertResult {
    assert(region.strand.has_value() && "Strand must be specified for insertion");

    auto &tree = featureTreeMap[region.referenceID];

    std::vector<size_t> indices;
    // Overlap with grace distance and blunt ends (+/- 1)
    tree.overlap(region.startPosition - graceDistance - 1, region.endPosition + graceDistance + 1,
                 indices);

    indices.erase(std::remove_if(indices.begin(), indices.end(),
                                 [&region, &tree](size_t index) {
                                     return tree.data(index).strand != *region.strand;
                                 }),
                  indices.end());

    if (indices.empty()) {
        return {insert(region), {}};
    }

    auto minStartIndex =
        *std::min_element(indices.begin(), indices.end(), [&tree](size_t lhs, size_t rhs) {
            return tree.data(lhs).startPosition < tree.data(rhs).startPosition;
        });

    auto maxEndIndex =
        *std::max_element(indices.begin(), indices.end(), [&tree](size_t lhs, size_t rhs) {
            return tree.data(lhs).endPosition < tree.data(rhs).endPosition;
        });

    dataTypes::Feature &minStartFeature = tree.data(minStartIndex);
    minStartFeature.startPosition = std::min(minStartFeature.startPosition, region.startPosition);
    minStartFeature.endPosition = std::max(tree.data(maxEndIndex).endPosition, region.endPosition);
    tree.setStart(minStartIndex, minStartFeature.startPosition);
    tree.setEnd(minStartIndex, minStartFeature.endPosition);

    std::vector<std::string> mergedFeatureIDs;
    mergedFeatureIDs.reserve(indices.size() - 1);

    for (unsigned long &indice : std::ranges::reverse_view(indices)) {
        if (indice != minStartIndex) {
            mergedFeatureIDs.push_back(std::move(tree.data(indice).id));
            tree.remove(indice);
        }
    }

    tree.index();
    return {.featureID = minStartFeature.id, .mergedFeatureIDs = std::move(mergedFeatureIDs)};
}

auto FeatureAnnotator::overlappingFeatures(const dataTypes::GenomicRegion &region,
                                           const Orientation orientation)
    -> std::vector<dataTypes::Feature> {
    std::vector<dataTypes::Feature> features;
    auto iterator = featureTreeMap.find(region.referenceID);
    if (iterator != featureTreeMap.end()) {
        std::vector<size_t> indices;
        iterator->second.overlap(region.startPosition, region.endPosition, indices);

        features.reserve(indices.size());

        for (const auto &index : indices) {
            const auto &feature = iterator->second.data(index);
            if ((orientation == Orientation::BOTH) ||
                (orientation == Orientation::OPPOSITE && feature.strand != *region.strand) ||
                (orientation == Orientation::SAME && feature.strand == *region.strand) ||
                (region.strand == std::nullopt)) {
                features.push_back(feature);
            }
        }
    }

    return features;
}

auto FeatureAnnotator::overlappingFeatureIterator(const dataTypes::GenomicRegion &region,
                                                  const Orientation orientation) const
    -> FeatureAnnotator::Results {
    std::vector<size_t> indices;

    auto iterator = featureTreeMap.find(region.referenceID);

    if (iterator == featureTreeMap.end()) [[unlikely]] {
        return {&iterator->second, indices, std::nullopt};
    }

    iterator->second.overlap(region.startPosition, region.endPosition, indices);

    std::optional<dataTypes::Strand> strand = std::nullopt;

    if (orientation == Orientation::SAME) {
        strand = region.strand;
    } else if (orientation == Orientation::OPPOSITE && region.strand.has_value()) {
        strand = !*region.strand;
    }

    return {&iterator->second, indices, strand};
}

auto FeatureAnnotator::getBestOverlappingFeature(const dataTypes::GenomicRegion &region,
                                                 const Orientation orientation) const
    -> std::optional<dataTypes::Feature> {
    auto overlapSizeWithRegion = [region](const dataTypes::Feature &feature) -> size_t {
        return std::min(region.endPosition, feature.endPosition) -
               std::max(region.startPosition, feature.startPosition);
    };

    auto featureIterator = overlappingFeatureIterator(region, orientation);

    auto maxOverlapSizeElement = std::max_element(
        featureIterator.begin(), featureIterator.end(),
        [&overlapSizeWithRegion](const dataTypes::Feature &lhs, const dataTypes::Feature &rhs) {
            return std::invoke(overlapSizeWithRegion, lhs) <
                   std::invoke(overlapSizeWithRegion, rhs);
        });

    if (maxOverlapSizeElement != featureIterator.end()) {
        return *maxOverlapSizeElement;
    }
    return std::nullopt;
}

auto FeatureAnnotator::getFeatureTreeMap() const -> const FeatureTreeMap & {
    return featureTreeMap;
}

auto FeatureAnnotator::getBestOverlappingFeature(
    const SamRecord &record, const std::deque<std::string> &referenceIDs,
    const Orientation orientation) const -> std::optional<dataTypes::Feature> {
    auto region = dataTypes::GenomicRegion::fromSamRecord(record, referenceIDs);

    if (!region.has_value()) {
        return std::nullopt;
    }

    return getBestOverlappingFeature(region.value(), orientation);
}

// Results and Iterator implementation
FeatureAnnotator::Results::Results(const IITree<int, dataTypes::Feature> *tree,
                                   const std::vector<size_t> &indices,
                                   const std::optional<dataTypes::Strand> strand)
    : tree(tree), indices(indices), strand(strand) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::begin() const {
    size_t startIndex = 0;
    if (strand.has_value()) {
        // Find the first index with the specified strand
        auto iterator = std::find_if(indices.begin(), indices.end(), [&](size_t index) {
            return tree->data(index).strand == *strand;
        });
        startIndex =
            iterator != indices.end() ? std::distance(indices.begin(), iterator) : indices.size();
    }
    return Iterator(tree, indices, startIndex, strand);
}

[[nodiscard]] auto FeatureAnnotator::Results::end() const -> FeatureAnnotator::Results::Iterator {
    return Iterator(tree, indices, indices.size(), strand);
}

FeatureAnnotator::Results::Iterator::Iterator(const IITree<int, dataTypes::Feature> *tree,
                                              const std::vector<size_t> &indices, size_t index,
                                              const std::optional<dataTypes::Strand> &strand)
    : tree(tree), indices(indices), current_index(index), strand(strand) {}

[[nodiscard]] auto FeatureAnnotator::Results::Iterator::operator*() const
    -> FeatureAnnotator::Results::Iterator::reference {
    if (current_index >= indices.size()) {
        throw std::out_of_range("Iterator out of range");
    }
    return tree->data(indices[current_index]);
}

[[nodiscard]] auto FeatureAnnotator::Results::Iterator::operator->() const
    -> FeatureAnnotator::Results::Iterator::pointer {
    if (current_index >= indices.size()) {
        throw std::out_of_range("Iterator out of range");
    }
    return &tree->data(indices[current_index]);
}

auto FeatureAnnotator::Results::Iterator::operator++() -> FeatureAnnotator::Results::Iterator & {
    auto matchStrand = [&](size_t index) {
        if (strand.has_value()) {
            return tree->data(indices[index]).strand == *strand;
        }
        return true;
    };

    ++current_index;

    while (current_index < indices.size() && !matchStrand(current_index)) {
        ++current_index;
    }

    return *this;
}

auto FeatureAnnotator::Results::Iterator::operator++(int) -> FeatureAnnotator::Results::Iterator {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}

}  // namespace annotation
