#include "Analyze.hpp"

#include <sys/stat.h>

#include <future>
#include <string>

#include "Constants.hpp"
#include "FeatureWriter.hpp"
#include "InteractionClusterGenerator.hpp"
#include "SplitRecordsParser.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {

// Analyze
void Analyze::process(const AnalyzeData &data) {
    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    std::vector<std::future<void>> futures;

    futures.reserve(data.treatmentSamples.size());
    for (const auto &sample : data.treatmentSamples) {
        futures.push_back(std::async(std::launch::async, &Analyze::processSample, this, sample));
    }

    if (!data.controlSamples.has_value()) {
        Logger::log(LogLevel::INFO, "No control samples provided");
        return;
    }

    for (const auto &sample : data.controlSamples.value()) {
        futures.push_back(std::async(std::launch::async, &Analyze::processSample, this, sample));
    }

    for (auto &fut : futures) {
        fut.get();
    }
}

void Analyze::processSample(AnalyzeSample sample) {
    Logger::log(LogLevel::INFO, "Processing sample: ", sample.input.sampleName);

    std::vector<InteractionCluster> clusters =
        SplitRecordsParser::parse(sample.input.splitAlignmentsPath);

    InteractionClusterGenerator clusterGenerator{sample.input.sampleName,
                                                 parameters.minimumClusterReadCount,
                                                 parameters.clusterDistanceThreshold};

    auto mergedClusters = clusterGenerator.mergeClusters(clusters);

    seqan3::sam_file_input splitsIn{sample.input.splitAlignmentsPath, sam_field_ids{}};

    auto &header = splitsIn.header();
    const std::deque<std::string> &referenceIDs = header.ref_ids();

    assignTranscriptsToClusters(mergedClusters, referenceIDs, sample);

    writeInteractionsToFile(mergedClusters, sample.input.sampleName, sample.output.interactionsPath,
                            referenceIDs);
}

void Analyze::assignTranscriptsToClusters(std::vector<InteractionCluster> &clusters,
                                          const std::deque<std::string> &referenceIDs,
                                          const AnalyzeSample &sample) {
    std::unordered_map<std::string, size_t> transcriptCounts;
    std::unordered_map<std::string, std::string> mergedTranscriptIDsIntoTranscriptIDs;

    annotation::FeatureAnnotator supplementaryFeatureAnnotator;

    const int mergeGraceDistance = parameters.clusterDistanceThreshold;
    const auto annotationOrientation = parameters.featureOrientation;

    transcriptCounts.reserve(clusters.size());

    auto processSegment = [&](const Segment &segment,
                              const InteractionCluster &cluster) -> std::string {
        const auto region = segment.toGenomicRegion(referenceIDs);

        auto updateTranscriptCounts =
            [&](const std::optional<dataTypes::Feature> &feature) -> std::optional<std::string> {
            if (feature.has_value()) {
                transcriptCounts[feature->id] += cluster.count;
                return feature->id;
            }
            return std::nullopt;
        };

        if (auto segmentFeature = updateTranscriptCounts(
                featureAnnotator.getBestOverlappingFeature(region, annotationOrientation))) {
            return *segmentFeature;
        }

        if (auto supplementarySegmentFeature =
                updateTranscriptCounts(supplementaryFeatureAnnotator.getBestOverlappingFeature(
                    region, annotationOrientation))) {
            return *supplementarySegmentFeature;
        }

        const auto mergeResult =
            supplementaryFeatureAnnotator.mergeInsert(region, mergeGraceDistance);
        transcriptCounts[mergeResult.featureID] += cluster.count;

        for (const auto &mergedFeatureID : mergeResult.mergedFeatureIDs) {
            mergedTranscriptIDsIntoTranscriptIDs[mergedFeatureID] = mergeResult.featureID;
            transcriptCounts[mergeResult.featureID] += transcriptCounts[mergedFeatureID];
            transcriptCounts.erase(mergedFeatureID);
        }

        return mergeResult.featureID;
    };

    for (auto &cluster : clusters) {
        const std::string firstSegmentTranscriptID =
            processSegment(cluster.segments.first, cluster);
        const std::string secondSegmentTranscriptID =
            processSegment(cluster.segments.second, cluster);

        cluster.transcriptIDs = std::make_pair(firstSegmentTranscriptID, secondSegmentTranscriptID);
    }

    for (auto &cluster : clusters) {
        if (cluster.transcriptIDs.has_value()) {
            auto &[firstTranscriptID, secondTranscriptID] = *cluster.transcriptIDs;

            if (auto iterator = mergedTranscriptIDsIntoTranscriptIDs.find(firstTranscriptID);
                iterator != mergedTranscriptIDsIntoTranscriptIDs.end()) {
                firstTranscriptID = iterator->second;
            }

            if (auto iterator = mergedTranscriptIDsIntoTranscriptIDs.find(secondTranscriptID);
                iterator != mergedTranscriptIDsIntoTranscriptIDs.end()) {
                secondTranscriptID = iterator->second;
            }
        }
    }

    annotation::FeatureWriter::write(supplementaryFeatureAnnotator.getFeatureTreeMap(),
                                     sample.output.supplementaryFeaturesPath,
                                     annotation::FileType::GFF);

    assignNonAnnotatedContiguousToSupplementaryFeatures(
        sample.input.unassignedContiguousAlignmentsPath, supplementaryFeatureAnnotator,
        transcriptCounts);

    assignAnnotatedContiguousFragmentCountsToTranscripts(
        sample.input.contiguousAlignmentsTranscriptCountsPath, transcriptCounts);

    const size_t totalFragmentCount =
        parseSampleFragmentCount(sample.input.sampleFragmentCountsPath);
    assignPValuesToClusters(clusters, transcriptCounts, totalFragmentCount);

    std::ofstream transcriptCountsOut(sample.output.interactionsTranscriptCountsPath);

    if (!transcriptCountsOut.is_open()) {
        Logger::log(LogLevel::ERROR,
                    "Could not open file: ", sample.output.interactionsTranscriptCountsPath);
    }

    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptCountsOut << transcriptID << "\t" << count << "\n";
    }
}

void Analyze::assignAnnotatedContiguousFragmentCountsToTranscripts(
    const fs::path &contiguousTranscriptCountsInPath,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    std::ifstream transcriptCountsIn(contiguousTranscriptCountsInPath);

    if (!transcriptCountsIn.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", contiguousTranscriptCountsInPath);
    }

    std::string line;
    while (std::getline(transcriptCountsIn, line)) {
        std::istringstream iss(line);

        size_t column = 0;
        std::string transcriptID;
        for (std::string token; std::getline(iss, token, '\t');) {
            if (column == 0) {
                transcriptID = token;
            } else if (column == 1) {
                transcriptCounts[transcriptID] += std::stoul(token);
            }
            ++column;
        }

        break;
    }
}

auto Analyze::getTranscriptProbabilities(
    const std::unordered_map<std::string, size_t> &transcriptCounts,
    const size_t totalTranscriptCount) -> std::unordered_map<std::string, double> {
    std::unordered_map<std::string, double> transcriptProbabilities;

    double cumulativeProbability = 0;
    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptProbabilities[transcriptID] =
            static_cast<double>(count) / static_cast<double>(totalTranscriptCount);
        cumulativeProbability += transcriptProbabilities[transcriptID];
    }

    // Normalize probabilities to sum to 1
    for (auto &[transcriptID, probability] : transcriptProbabilities) {
        probability /= cumulativeProbability;
    }

    return transcriptProbabilities;
}

void Analyze::assignPAdjustedValuesToClusters(std::vector<InteractionCluster> &clusters) {
    std::vector<std::pair<double, size_t>> pValuesWithIndex;
    pValuesWithIndex.reserve(clusters.size());  // Reserve capacity to avoid reallocations

    for (size_t i = 0; i < clusters.size(); ++i) {
        if (clusters[i].pValue.has_value()) {
            pValuesWithIndex.emplace_back(clusters[i].pValue.value(), i);
        }
    }

    std::sort(pValuesWithIndex.begin(), pValuesWithIndex.end());

    const auto numPValues = static_cast<double>(pValuesWithIndex.size());
    double minAdjPValue = 1.0;
    for (size_t i = pValuesWithIndex.size(); i-- > 0;) {
        auto rank = static_cast<double>(i + 1);
        double adjustedPValue =
            std::min(minAdjPValue, (pValuesWithIndex[i].first * numPValues / rank));
        minAdjPValue = adjustedPValue;
        clusters[pValuesWithIndex[i].second].pAdj = adjustedPValue;
    }
}

void Analyze::assignPValuesToClusters(
    std::vector<InteractionCluster> &clusters,
    const std::unordered_map<std::string, size_t> &transcriptCounts,
    const size_t totalTranscriptCount) {
    const auto transcriptProbabilities =
        getTranscriptProbabilities(transcriptCounts, totalTranscriptCount);

    for (auto &cluster : clusters) {
        if (!cluster.transcriptIDs.has_value()) {
            Logger::log(LogLevel::WARNING, "Some cluster is missing transcript IDs");
            continue;
        }

        const auto &firstTranscriptID = cluster.transcriptIDs->first;
        const auto &secondTranscriptID = cluster.transcriptIDs->second;

        auto findProbability = [&](const std::string &transcriptID) -> std::optional<double> {
            auto iterator = transcriptProbabilities.find(transcriptID);
            if (iterator != transcriptProbabilities.end()) {
                return iterator->second;
            }
            return std::nullopt;
        };

        auto firstTranscriptProbability = findProbability(firstTranscriptID);
        auto secondTranscriptProbability = findProbability(secondTranscriptID);

        if (!firstTranscriptProbability || !secondTranscriptProbability) {
            Logger::log(LogLevel::WARNING, "Could not find transcript probabilities for cluster.");
            continue;
        }

        const double ligationByChanceProbability =
            (firstTranscriptID == secondTranscriptID)
                ? firstTranscriptProbability.value() * secondTranscriptProbability.value()
                : 2 * firstTranscriptProbability.value() * secondTranscriptProbability.value();
        ;

        const auto binomialDistribution =
            math::binomial_distribution((double)totalTranscriptCount, ligationByChanceProbability);

        const double pValue = 1 - math::cdf(binomialDistribution, cluster.count);

        cluster.pValue = pValue;
    }

    assignPAdjustedValuesToClusters(clusters);
}

void Analyze::assignNonAnnotatedContiguousToSupplementaryFeatures(
    const fs::path &unassignedSingletonsInPath, annotation::FeatureAnnotator &featureAnnotator,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    seqan3::sam_file_input unassignedSingletonsIn{unassignedSingletonsInPath.string(),
                                                  sam_field_ids{}};

    const auto annotationOrientation = parameters.featureOrientation;

    for (auto &&records : unassignedSingletonsIn) {
        const auto region = dataTypes::GenomicRegion::fromSamRecord(
            records, unassignedSingletonsIn.header().ref_ids());

        if (!region) {
            continue;
        }

        const auto bestFeature =
            featureAnnotator.getBestOverlappingFeature(region.value(), annotationOrientation);

        if (!bestFeature) {
            continue;
        }

        const std::string &transcriptID = bestFeature->id;
        assert(transcriptCounts.contains(transcriptID));
        ++transcriptCounts[transcriptID];
    }
}

void Analyze::writeBEDLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                                 const std::deque<std::string> &referenceIDs,
                                 const std::string &color, std::ofstream &bedOut) {
    bedOut << getReferenceID(cluster.segments.first.referenceIDIndex, referenceIDs) << "\t";
    bedOut << cluster.segments.first.start << "\t";
    bedOut << cluster.segments.first.end << "\t";
    bedOut << "cluster" << clusterID << "_segment1"
           << "\t";
    bedOut << "0\t";
    bedOut << static_cast<char>(cluster.segments.first.strand) << "\t";
    bedOut << cluster.segments.first.start << "\t";
    bedOut << cluster.segments.first.end << "\t";
    bedOut << color << "\n";

    bedOut << getReferenceID(cluster.segments.second.referenceIDIndex, referenceIDs) << "\t";
    bedOut << cluster.segments.second.start << "\t";
    bedOut << cluster.segments.second.end << "\t";
    bedOut << "cluster" << clusterID << "_segment2"
           << "\t";
    bedOut << "0\t";
    bedOut << static_cast<char>(cluster.segments.second.strand) << "\t";
    bedOut << cluster.segments.second.start << "\t";
    bedOut << cluster.segments.second.end << "\t";
    bedOut << color << "\n";
}

void Analyze::writeBEDArcLineToFile(const InteractionCluster &cluster, const std::string &clusterID,
                                    const std::deque<std::string> &referenceIDs,
                                    std::ofstream &bedOut) {
    bedOut << getReferenceID(cluster.segments.first.referenceIDIndex, referenceIDs) << "\t";
    bedOut << cluster.segments.first.start << "\t";
    bedOut << cluster.segments.second.start << "\t";
    bedOut << "cluster" << clusterID << "\n";
}

void Analyze::writeInteractionLineToFile(const InteractionCluster &cluster,
                                         const std::string &clusterID,
                                         const std::deque<std::string> &referenceIDs,
                                         std::ofstream &interactionOut) {
    interactionOut << "cluster" << clusterID << "\t";

    interactionOut << cluster.transcriptIDs->first << "\t";
    interactionOut << getReferenceID(cluster.segments.first.referenceIDIndex, referenceIDs) << "\t";
    interactionOut << static_cast<char>(cluster.segments.first.strand) << "\t";
    interactionOut << cluster.segments.first.start << "\t";
    interactionOut << cluster.segments.first.end << "\t";

    interactionOut << cluster.transcriptIDs->second << "\t";
    interactionOut << getReferenceID(cluster.segments.second.referenceIDIndex, referenceIDs)
                   << "\t";
    interactionOut << static_cast<char>(cluster.segments.second.strand) << "\t";
    interactionOut << cluster.segments.second.start << "\t";
    interactionOut << cluster.segments.second.end << "\t";

    interactionOut << cluster.count << "\t";
    interactionOut << cluster.complementarityStatistics() << "\t";
    interactionOut << cluster.hybridizationEnergyStatistics() << "\t";
    interactionOut << cluster.pValue.value_or(1.0) << "\t";
    interactionOut << cluster.pAdj.value_or(1.0);
    interactionOut << "\n";
}

void Analyze::writeInteractionsToFile(const std::vector<InteractionCluster> &mergedClusters,
                                      const std::string &sampleName, const fs::path &clusterOutPath,
                                      const std::deque<std::string> &referenceIDs) const {
    std::ofstream interactionOut(clusterOutPath.string());
    if (!interactionOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", clusterOutPath.string());
    }

    interactionOut << "cluster_ID\tfst_feat_id\tfst_seg_chr\tfst_seg_strd\tfst_seg_strt\tfst_seg_"
                      "end\tsec_feat_id\t"
                      "sec_seg_chr\tsec_seg_strd\tsec_seg_strt\tsec_seg_end\tno_splits\t"
                      "gcs\tghs\tp_value\tpadj_value\n";

    fs::path bedOutPath = clusterOutPath;
    bedOutPath.replace_extension(".bed");

    std::ofstream bedOut(bedOutPath.string());
    if (!bedOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", bedOutPath.string());
    }

    bedOut << "track name=\"" << sampleName
           << " RNA-RNA interactions\" description=\"Segments of interacting RNA "
              "clusters derived from DDD-Experiment\" itemRgb=\"On\"\n";

    fs::path arcOutBase = clusterOutPath.stem();
    arcOutBase += "_arc";
    fs::path arcOutPath = clusterOutPath.parent_path() / arcOutBase;
    arcOutPath.replace_extension(".bed");

    std::ofstream arcOut(arcOutPath.string());
    if (!arcOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", arcOutPath.string());
    }

    arcOut << "track graphType=arc name=\"" << sampleName
           << " RNA-RNA interactions arcs\" description=\"Segments of interacting RNA "
              "clusters derived from DDD-Experiment\" itemRgb=\"On\"\n";

    const double pValueThreshold = parameters.padjThreshold;
    const auto clusterCountThreshold = parameters.minimumClusterReadCount;

    auto isInteractionValid = [&](const InteractionCluster &cluster) {
        return cluster.pAdj.has_value() && cluster.pAdj.value() <= pValueThreshold &&
               cluster.count >= static_cast<int>(clusterCountThreshold);
    };

    size_t intramolecularCount = 0;
    size_t intermolecularCount = 0;

    size_t clusterID = 0;
    for (const auto &interaction : mergedClusters) {
        if (!isInteractionValid(interaction)) {
            continue;
        }

        if (interaction.transcriptIDs.has_value(),
            interaction.transcriptIDs.value().first == interaction.transcriptIDs.value().second) {
            ++intramolecularCount;
        } else {
            ++intermolecularCount;
        }

        const std::string randomColor = helper::generateRandomHexColor();
        writeBEDLineToFile(interaction, std::to_string(clusterID), referenceIDs, randomColor,
                           bedOut);
        writeBEDArcLineToFile(interaction, std::to_string(clusterID), referenceIDs, arcOut);
        writeInteractionLineToFile(interaction, std::to_string(clusterID), referenceIDs,
                                   interactionOut);
        ++clusterID;
    }
    Logger::log(LogLevel::INFO, "(", sampleName, ") After filtering kept ",
                intramolecularCount + intermolecularCount, " split interactions");
    Logger::log(LogLevel::INFO, "(", sampleName, ") Of which ", intramolecularCount,
                " are intramolecular interactions");
    Logger::log(LogLevel::INFO, "(", sampleName, ") Of which ", intermolecularCount,
                " are intermolecular interactions");
}

auto Analyze::parseSampleFragmentCount(const fs::path &sampleCountsInPath) -> size_t {
    std::ifstream sampleCountsIn(sampleCountsInPath);

    if (!sampleCountsIn.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", sampleCountsInPath.string());
    }

    sampleCountsIn.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip header

    size_t totalTranscriptCount = 0;
    std::string line;
    while (std::getline(sampleCountsIn, line)) {
        std::istringstream iss(line);

        size_t column = 0;
        for (std::string token; std::getline(iss, token, '\t');) {
            if (column != 0) {
                totalTranscriptCount += std::stoul(token);
            }
            ++column;
        }

        break;  // Should only be one line
    }

    return totalTranscriptCount;
}

}  // namespace pipelines::analyze
