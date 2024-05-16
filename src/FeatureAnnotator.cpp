#include "FeatureAnnotator.hpp"

using namespace Annotation;

FeatureAnnotator::FeatureAnnotator(fs::path featureFilePath,
                                   const std::unordered_set<std::string> &includedFeatures,
                                   const std::optional<std::string> &featureIDFlag = std::nullopt) {
    dtp::FeatureMap featureMap =
        FeatureParser(includedFeatures, featureIDFlag).parse(featureFilePath);

    for (const auto &[seqid, features] : featureMap) {
        IITree<int, dtp::Feature> tree;
        for (const auto &feature : features) {
            tree.add(feature.start, feature.end, feature);
        }
        featureTreeMap[seqid] = tree;
        featureTreeMap[seqid].index();
    }
}

std::vector<dtp::Feature> FeatureAnnotator::overlappingFeatures(
    const dtp::GenomicRegion &region) const {
    std::vector<dtp::Feature> features;
    if (featureTreeMap.contains(region.seqid)) {
        std::vector<size_t> indices;
        featureTreeMap.at(region.seqid).overlap(region.start, region.end, indices);
        for (const auto &index : indices) {
            features.push_back(featureTreeMap.at(region.seqid).data(index));
        }
    }
    return features;
}

FeatureAnnotator::Results FeatureAnnotator::getOverlappingFeatureIterator(
    const dtp::GenomicRegion &region) {
    if (auto it = featureTreeMap.find(region.seqid); it != featureTreeMap.end()) {
        std::vector<size_t> indices;
        it->second.overlap(region.start, region.end, indices);
        return Results(it->second, std::move(indices));
    }
    return Results(IITree<int, dtp::Feature>{}, {});
}

FeatureAnnotator::Results::Results(const IITree<int, dtp::Feature> &tree,
                                   const std::vector<size_t> indices)
    : tree_(tree), indices_(indices) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::begin() const {
    return Iterator(tree_, indices_, 0);
}

[[nodiscard]] FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::end() const {
    return Iterator(tree_, indices_, indices_.size());
}

FeatureAnnotator::Results::Iterator::Iterator(const IITree<int, dtp::Feature> &tree,
                                              const std::vector<size_t> indices, size_t index)
    : tree_(tree), indices_(indices), current_index_(index) {}

[[nodiscard]] FeatureAnnotator::Results::Iterator::reference
FeatureAnnotator::Results::Iterator::operator*() const {
    return tree_.data(indices_[current_index_]);
}

[[nodiscard]] FeatureAnnotator::Results::Iterator::pointer
FeatureAnnotator::Results::Iterator::operator->() const {
    return &tree_.data(indices_[current_index_]);
}

FeatureAnnotator::Results::Iterator &FeatureAnnotator::Results::Iterator::operator++() {
    ++current_index_;
    return *this;
}

FeatureAnnotator::Results::Iterator FeatureAnnotator::Results::Iterator::operator++(int) {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}