#include "Detect.hpp"

Detect::Detect(po::variables_map params)
    : params(params),
      minReadLength(params["minlen"].as<std::size_t>()),
      minMapQuality(params["mapqmin"].as<int>()),
      excludeSoftClipping(params["exclclipping"].as<bool>()),
      annotationOrientation(params["orientation"].as<Annotation::Orientation>()),
      featureAnnotator(params["features"].as<std::string>(),
                       params["featuretypes"].as<std::vector<std::string>>()),
      splitRecordsEvaluator(SplitRecordsEvaluator(getSplitRecordsEvaluatorParameters(params))) {}

const std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                   SplitRecordsEvaluationParameters::SplicingParameters>
Detect::getSplitRecordsEvaluatorParameters(const po::variables_map &params) const {
    if (params["splicing"].as<bool>()) {
        return SplitRecordsEvaluationParameters::SplicingParameters{
            .baseParameters = {.minComplementarity = params["cmplmin"].as<double>(),
                               .minComplementarityFraction = params["sitelenratio"].as<double>(),
                               .mfeThreshold = params["mfe"].as<double>()},
            .orientation = params["orientation"].as<Annotation::Orientation>(),
            .splicingTolerance = params["splicingtolerance"].as<int>(),
            .featureAnnotator = featureAnnotator};
    } else {
        return SplitRecordsEvaluationParameters::BaseParameters{
            .minComplementarity = params["cmplmin"].as<double>(),
            .minComplementarityFraction = params["sitelenratio"].as<double>(),
            .mfeThreshold = params["mfe"].as<double>()};
    }
}

std::deque<std::string> const &Detect::getReferenceIDs(const fs::path &mappingsInPath) const {
    seqan3::sam_file_input alignmentsIn{mappingsInPath.string(), sam_field_ids{}};
    return alignmentsIn.header().ref_ids();
}

Detect::Results Detect::iterateSortedMappingsFile(
    const std::string &mappingsInPath, const std::string &splitsPath,
    const std::string &multSplitsPath, const fs::path &unassignedSingletonRecordsOutPath) {
    seqan3::sam_file_input alignmentsIn{mappingsInPath, sam_field_ids{}};

    auto &header = alignmentsIn.header();
    std::vector<size_t> referenceLengths{};
    std::ranges::transform(header.ref_id_info, std::back_inserter(referenceLengths),
                           [](auto const &info) { return std::get<0>(info); });
    const std::deque<std::string> &refIDs = header.ref_ids();

    seqan3::sam_file_output splitsOut{splitsPath, refIDs, referenceLengths, sam_field_ids{}};
    seqan3::sam_file_output multiSplitsOut{multSplitsPath, refIDs, referenceLengths,
                                           sam_field_ids{}};

    seqan3::sam_file_output unassignedSingletonRecordsOut{
        unassignedSingletonRecordsOutPath.string(), refIDs, referenceLengths};

    std::string currentRecordID;
    std::vector<SamRecord> currentRecordList;
    std::unordered_set<size_t> invalidRecordHitGroups;

    size_t recordsCount = 0;
    size_t totalSplitFragmentsCount = 0;
    size_t totalSingletonFragmentsCount = 0;

    TranscriptCounts singletonTranscriptCounts;

    auto assignSingletonTranscriptCount = [&](SamRecord &record) {
        if (!record.reference_id() || !record.reference_position()) [[unlikely]] {
            return;
        }

        const auto region = GenomicRegion::fromSamRecord(record, refIDs);

        if (!region) {
            return;
        }

        const auto bestFeature =
            featureAnnotator.getBestOverlappingFeature(region.value(), annotationOrientation);

        if (bestFeature) {
            singletonTranscriptCounts[bestFeature->id]++;
        } else {
            unassignedSingletonRecordsOut.push_back(record);
        }
    };

    for (auto &record : alignmentsIn) {
        recordsCount++;

        // This depends on the SAM file being sorted by read ID
        bool isNewRecordID = !currentRecordID.empty() && (currentRecordID != record.id());

        if (isNewRecordID) {
            // Remove current records that belong to invalid hit groups
            std::erase_if(currentRecordList, [&](const SamRecord &record) {
                return invalidRecordHitGroups.contains(record.tags().get<"HI"_tag>());
            });

            size_t fragmentsCount =
                processReadRecords(currentRecordList, refIDs, splitsOut, multiSplitsOut);
            totalSplitFragmentsCount += fragmentsCount;

            currentRecordList.clear();
            invalidRecordHitGroups.clear();
        }

        currentRecordID = record.id();

        // Read is completely invalid and does not meet minimum requirements -> Not counted
        if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped) ||
            record.mapping_quality() < minMapQuality || record.sequence().size() < minReadLength) {
            invalidRecordHitGroups.insert(record.tags().get<"HI"_tag>());
            continue;
        }

        // Read is singleton and does not contain any split information -> Counted for total
        // fragments
        if (!record.tags().contains("XJ"_tag) || record.tags().get<"XJ"_tag>() < 2) {
            if (static_cast<bool>(record.flag() & seqan3::sam_flag::secondary_alignment)) {
                continue;
            }

            assignSingletonTranscriptCount(record);
            totalSingletonFragmentsCount++;

            continue;
        }

        // Read is part of a split read and valid but a previous record with the same hit group
        // was invalid -> Skip
        if (invalidRecordHitGroups.contains(record.tags().get<"HI"_tag>())) {
            continue;
        }

        currentRecordList.push_back(std::move(record));
    }

    size_t fragmentsCount =
        processReadRecords(currentRecordList, refIDs, splitsOut, multiSplitsOut);
    totalSplitFragmentsCount += fragmentsCount;

    Logger::log(LogLevel::INFO, "Processed ", recordsCount, " reads. Found ",
                totalSplitFragmentsCount, " split fragments and ", totalSingletonFragmentsCount,
                " singleton fragments.");

    return {singletonTranscriptCounts, totalSplitFragmentsCount, totalSingletonFragmentsCount};
}

/**
 * Retrieves the final optimal split records and returns the count of fragments for a specific
 * record id. Expects read records to stem exactly from one record id.
 *
 * @param readRecords The vector of SamRecord objects representing the read records.
 * @param splitsOut The output stream for split records.
 * @param multiSplitsOut The output stream for multi split records.
 * @return The count of fragments for a specific record id.
 */
size_t Detect::processReadRecords(const std::vector<SamRecord> &readRecords,
                                  const std::deque<std::string> &referenceIDs, auto &splitsOut,
                                  auto &multiSplitsOut) {
    if (readRecords.empty()) {
        return 0;
    }

    const auto splitRecords = getSplitRecords(readRecords, referenceIDs);

    if (!splitRecords.has_value()) {
        return 0;
    }

    // TODO Implement multi split writing
    writeSamFile(splitsOut, splitRecords.value().splitRecords);

    return splitRecords.value().splitRecords.size();
}

std::optional<SplitRecords> Detect::constructSplitRecords(const SamRecord &readRecord) {
    // Number of expected split records for within the record
    const size_t expectedSplitRecords = readRecord.tags().get<"XH"_tag>();

    SplitRecords splitRecords{};
    splitRecords.reserve(expectedSplitRecords);

    std::vector<seqan3::cigar> currentCigar{};
    splitRecords.reserve(10);

    size_t referencePosition = readRecord.reference_position().value_or(0);
    size_t startPosRead{};
    size_t endPosRead{};  // absolute position in read/alignment (e.g., 1 to end)
    size_t startPosSplit{};
    size_t endPosSplit{};  // position in split read (e.g., XX:i to XY:i / 14 to 20)

    bool isValid = true;
    int nextSplitReferenceShift = 0;

    auto const addSplitRecord = [&]() {
        const auto splitSeq =
            readRecord.sequence() | seqan3::views::slice(startPosRead, endPosRead);
        const auto splitQual =
            readRecord.base_qualities() | seqan3::views::slice(startPosRead, endPosRead);

        if (splitSeq.size() < minReadLength) {
            isValid = false;
            return;
        }

        seqan3::sam_tag_dictionary tags{};
        tags.get<"XX"_tag>() = startPosSplit;
        tags.get<"XY"_tag>() = endPosSplit;
        tags.get<"XN"_tag>() = splitRecords.size();

        splitRecords.emplace_back(readRecord.id(), readRecord.flag(), readRecord.reference_id(),
                                  referencePosition, readRecord.mapping_quality(), currentCigar,
                                  seqan3::dna5_vector(splitSeq.begin(), splitSeq.end()),
                                  std::vector<seqan3::phred42>(splitQual.begin(), splitQual.end()),
                                  std::move(tags));
    };

    auto const addOtherCigar = [&](const auto &cigar) {
        const auto cigarValue = get<0>(cigar);
        endPosRead += cigarValue;
        endPosSplit += cigarValue;
        currentCigar.push_back(cigar);
    };

    auto const addInsertionCigar = [&](const auto &cigar) {
        const auto cigarValue = get<0>(cigar);
        endPosRead += cigarValue;
        endPosSplit += cigarValue;
        nextSplitReferenceShift -= cigarValue;
        currentCigar.push_back(cigar);
    };

    auto const addDeletionCigar = [&](const auto &cigar) {
        currentCigar.push_back(cigar);
        nextSplitReferenceShift += get<0>(cigar);
    };

    auto const addSoftClipCigar = [&](const auto &cigar) {
        if (!excludeSoftClipping) return addOtherCigar(cigar);

        /* If current cigar is empty, we are at the beginning of the read in case of soft
        clipping at the end of the read it is just ignored */
        if (currentCigar.empty()) {
            const auto cigarValue = get<0>(cigar);
            nextSplitReferenceShift -= cigarValue;
            startPosRead += cigarValue;
            endPosRead += cigarValue;
            startPosSplit += cigarValue;
            endPosSplit += cigarValue;
        }
    };

    auto const addSkipCigar = [&](const auto &cigar) {
        if (currentCigar.empty()) {
            return;
        }

        addSplitRecord();

        // Set up positions for the next split
        const auto cigarValue = get<0>(cigar);
        assert((cigarValue + endPosRead + nextSplitReferenceShift) >= 0);
        referencePosition += cigarValue + endPosRead + nextSplitReferenceShift;
        startPosSplit = endPosSplit;
        startPosRead = endPosRead;
        nextSplitReferenceShift = 0;
        currentCigar.clear();
    };

    for (const auto &cigar : readRecord.cigar_sequence()) {
        if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
            cigar == 'X'_cigar_operation) {
            addOtherCigar(cigar);
        } else if (cigar == 'I'_cigar_operation) {
            addInsertionCigar(cigar);
        } else if (cigar == 'D'_cigar_operation) {
            addDeletionCigar(cigar);
        } else if (cigar == 'S'_cigar_operation) {
            addSoftClipCigar(cigar);
        } else if (cigar == 'N'_cigar_operation) {
            addSkipCigar(cigar);
        }
    }

    addSplitRecord();

    if (!isValid) {
        return std::nullopt;
    }

    if (splitRecords.size() != expectedSplitRecords) {
        Logger::log(LogLevel::WARNING, "Expected ", expectedSplitRecords,
                    " split records, but got ", splitRecords.size(),
                    ". Record ID: ", readRecord.id());

        return std::nullopt;
    }

    return splitRecords;
}

/**
 * Constructs split records for a list of records with one or more elements that comprise a
 * splitted read.
 *
 * @param readRecords The vector of read records.
 * @return An optional containing the split records if construction is successful, otherwise
 * std::nullopt.
 */
std::optional<SplitRecords> Detect::constructSplitRecords(
    const std::vector<SamRecord> &readRecords) {
    // Number of expected split records for whole read
    const size_t expectedSplitRecords = readRecords.front().tags().get<"XJ"_tag>();

    SplitRecords splitRecords{};
    splitRecords.reserve(expectedSplitRecords);

    for (const auto &record : readRecords) {
        auto splitRecord = constructSplitRecords(record);

        if (!splitRecord.has_value()) {
            return std::nullopt;
        }

        splitRecords.insert(splitRecords.end(), splitRecord->begin(), splitRecord->end());
    }

    if (splitRecords.size() != expectedSplitRecords) {
        Logger::log(LogLevel::WARNING, "Expected ", expectedSplitRecords,
                    " split records, but got ", splitRecords.size(),
                    ". Record ID: ", readRecords.front().id());

        return std::nullopt;
    }

    return splitRecords;
}

/**
 * Retrieves the split records from the given vector of read records.
 *
 * This function groups the read records based on their "HI" tag and constructs split records
 * for each group. It then evaluates the split records and returns the best evaluated split
 * records.
 *
 * @param readRecords The vector of read records.
 * @return An optional containing the best evaluated split records, or an empty optional if no
 *         split records were found.
 */
std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords> Detect::getSplitRecords(
    const std::vector<SamRecord> &readRecords, const std::deque<std::string> &referenceIDs) {
    std::unordered_map<size_t, std::vector<SamRecord>> recordHitGroups{};
    for (const auto &record : readRecords) {
        recordHitGroups[record.tags().get<"HI"_tag>()].push_back(record);
    }

    std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords> bestSplitRecords{};

    const auto insertBestSplitRecords = [&](SplitRecords &splitRecords) {
        const auto evaluationResult = splitRecordsEvaluator.evaluate(splitRecords, referenceIDs);

        if (std::holds_alternative<SplitRecordsEvaluator::EvaluatedSplitRecords>(
                evaluationResult)) {
            const auto &evaluatedSplitRecords =
                std::get<SplitRecordsEvaluator::EvaluatedSplitRecords>(evaluationResult);

            if (!bestSplitRecords.has_value() || evaluatedSplitRecords > bestSplitRecords.value()) {
                bestSplitRecords.emplace(
                    std::get<SplitRecordsEvaluator::EvaluatedSplitRecords>(evaluationResult));
            }
        } else {
            Logger::log(LogLevel::DEBUG, "Split records failed evaluation. Reason: ",
                        std::get<SplitRecordsEvaluator::FilterReason>(evaluationResult));
        }
    };

    for (const auto &[_, hitGroup] : recordHitGroups) {
        auto splitRecords = constructSplitRecords(hitGroup);

        if (splitRecords.has_value()) {
            insertBestSplitRecords(splitRecords.value());
        }
    }

    return bestSplitRecords;
}

void Detect::writeSamFile(auto &samOut, const std::vector<SamRecord> &splitRecords) {
    for (auto &&record : splitRecords) {
        auto [id, flag, ref_id, ref_offset, mapq, cigar, seq, qual, tags] = record;
        samOut.emplace_back(id, flag, ref_id, ref_offset, mapq, cigar, seq, qual, tags);
    }
}

std::vector<Detect::ChunkedInOutFilePaths> Detect::prepareInputOutputFiles(
    const fs::path &mappingsFilePath, const fs::path &splitsFilePath,
    const int mappingRecordsCount) {
    // create tmp folder (for inputs)
    fs::path tmpInDir = mappingsFilePath.parent_path() / "tmpIn";
    fs::path tmpSplitsOutDir = splitsFilePath.parent_path() / "tmpSplit";
    fs::path tmpMultiSplitsOutDir = splitsFilePath.parent_path() / "tmpMsplit";
    fs::path tmpUnassignedSingletonRecordsOutDirPath =
        splitsFilePath.parent_path() / "tmpUnassignedSingleton";

    helper::createTmpDir(tmpInDir);
    helper::createTmpDir(tmpSplitsOutDir);
    helper::createTmpDir(tmpMultiSplitsOutDir);
    helper::createTmpDir(tmpUnassignedSingletonRecordsOutDirPath);

    std::string tmpInPath = (tmpInDir / "tmp").string();

    const std::vector<fs::path> splitMappingsFilePaths =
        splitMappingsFile(mappingsFilePath, tmpInPath, mappingRecordsCount);

    std::vector<Detect::ChunkedInOutFilePaths> subFiles;
    for (auto &file : splitMappingsFilePaths) {
        fs::path tmpSplitsOutPath = tmpSplitsOutDir / file.filename();
        fs::path tmpMultiSplitsOutPath = tmpMultiSplitsOutDir / file.filename();
        fs::path tmpUnassignedSingletonRecordsOutPath =
            tmpUnassignedSingletonRecordsOutDirPath / file.filename();
        subFiles.push_back(
            {file, tmpSplitsOutPath, tmpMultiSplitsOutPath, tmpUnassignedSingletonRecordsOutPath});
    }

    return subFiles;
}

std::vector<fs::path> Detect::splitMappingsFile(const fs::path &mappingsFilePath,
                                                const fs::path &tmpInPath, const int entries) {
    std::ifstream inputFile(mappingsFilePath, std::ios::in);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open the input file.");
    }

    std::ofstream outputFile;
    std::string line, prev;
    std::stringstream hdrStream;
    int numLines = 0, numBlocks = 0;

    std::string outPathPrefix = tmpInPath.string() + "_";
    std::vector<fs::path> tmpOutFiles;

    while (std::getline(inputFile, line)) {
        if (line.starts_with('@')) {
            hdrStream << line << '\n';
            continue;
        }

        ++numLines;
        if ((numLines % entries) == 1) {
            if (line == prev) {
                --numLines;
            } else {
                if (outputFile.is_open()) {
                    outputFile.close();
                }
                std::string outFileName = outPathPrefix + std::to_string(++numBlocks) + ".sam";
                tmpOutFiles.push_back(outFileName);
                outputFile.open(outFileName, std::ios::out);
                if (!outputFile.is_open()) {
                    throw std::runtime_error("Could not open the output file.");
                }
                outputFile << hdrStream.str();
                numLines = 1;
            }
        }
        outputFile << line << '\n';
        prev = line;
    }

    if (outputFile.is_open()) {
        outputFile.close();
    }

    inputFile.close();

    return tmpOutFiles;
}

void Detect::writeStatisticsFile(const Results &results, const std::string &sampleName,
                                 const fs::path &statsFilePath) const {
    std::ofstream statsFileStream(statsFilePath, std::ios::app);

    if (!statsFileStream.is_open()) {
        throw std::runtime_error("Could not open the stats file.");
    }

    statsFileStream << "sample\tsplits\tsingletons" << std::endl;
    statsFileStream << sampleName << "\t" << results.splitFragmentsCount << "\t"
                    << results.singletonFragmentsCount << "\n";
}

void Detect::mergeResults(Detect::Results &results, const Detect::Results &newResults) const {
    results.splitFragmentsCount += newResults.splitFragmentsCount;
    results.singletonFragmentsCount += newResults.singletonFragmentsCount;

    for (const auto &[transcript, count] : newResults.transcriptCounts) {
        results.transcriptCounts[transcript] += count;
    }
}

void Detect::writeTranscriptCountsFile(const fs::path &transcriptCountsFilePath,
                                       const TranscriptCounts &transcriptCounts) const {
    std::ofstream transcriptCountsFileStream(transcriptCountsFilePath);

    if (!transcriptCountsFileStream.is_open()) {
        throw std::runtime_error("Could not open the transcript counts file.");
    }

    for (const auto &[transcript, count] : transcriptCounts) {
        transcriptCountsFileStream << transcript << "\t" << count << "\n";
    }
}

void Detect::start(pt::ptree sample) {
    pt::ptree input = sample.get_child("input");
    fs::path mappingRecordsInPath = input.get<std::string>("matched");

    int threads = params["threads"].as<int>();

    if (threads < 1) {
        throw std::runtime_error("Number of threads must be at least 1.");
    }

    size_t entries = std::round(helper::countSamEntries(mappingRecordsInPath) / threads);

    // output
    pt::ptree output = sample.get_child("output");
    fs::path mergedSplitsOutPath = output.get<std::string>("splits");
    fs::path mergedMultiSplitsOutPath = output.get<std::string>("multsplits");
    fs::path unassignedSingletonRecordsOutPath = output.get<std::string>("singletonUnassigned");
    fs::path singletonTranscriptCountsOutPath = output.get<std::string>("singleton");

    std::vector<ChunkedInOutFilePaths> subFiles =
        prepareInputOutputFiles(mappingRecordsInPath, mergedSplitsOutPath, entries);

    Detect::Results detectResults;

#pragma omp parallel for num_threads(threads)
    for (const auto &subFileChunk : subFiles) {
        const auto newResults = iterateSortedMappingsFile(
            subFileChunk.mappingsInPath.string(), subFileChunk.splitsOutPath.string(),
            subFileChunk.multSplitsOutPath.string(),
            subFileChunk.unassignedSingletonRecordsOutPath);

#pragma omp critical
        { mergeResults(detectResults, newResults); }
    }

    writeTranscriptCountsFile(singletonTranscriptCountsOutPath, detectResults.transcriptCounts);

    // Merge results from all processing chunks
    std::vector<fs::path> splitsChunkPaths;
    std::vector<fs::path> multiSplitsChunkPaths;
    std::vector<fs::path> unassignedSingletonRecordsChunkPaths;

    for (const auto &subFileChunk : subFiles) {
        splitsChunkPaths.push_back(subFileChunk.splitsOutPath);
        multiSplitsChunkPaths.push_back(subFileChunk.multSplitsOutPath);
        unassignedSingletonRecordsChunkPaths.push_back(
            subFileChunk.unassignedSingletonRecordsOutPath);
    }

    helper::mergeSamFiles(splitsChunkPaths, mergedSplitsOutPath);
    helper::mergeSamFiles(multiSplitsChunkPaths, mergedMultiSplitsOutPath);
    helper::mergeSamFiles(unassignedSingletonRecordsChunkPaths, unassignedSingletonRecordsOutPath);

    // Remove temporary files
    fs::remove_all(subFiles.front().mappingsInPath.parent_path());
    fs::remove_all(subFiles.front().splitsOutPath.parent_path());
    fs::remove_all(subFiles.front().multSplitsOutPath.parent_path());
    fs::remove_all(subFiles.front().unassignedSingletonRecordsOutPath.parent_path());

    fs::path statsFilePath = output.get<std::string>("stats");
    writeStatisticsFile(detectResults, mergedSplitsOutPath.stem().string(), statsFilePath);
}
