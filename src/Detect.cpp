#include "Detect.hpp"

using namespace seqan3::literals;

Detect::Detect(po::variables_map params)
    : params(params),
      minReadLength(params["minlen"].as<std::size_t>()),
      excludeSoftClipping(params["exclclipping"].as<bool>()),
      minComplementarity(params["cmplmin"].as<double>()),
      minFraction(params["sitelenratio"].as<double>()) {}
//
void Detect::iterateSortedMappingsFile(const std::string &mappingsInPath,
                                       const std::string &splitsPath,
                                       const std::string &multSplitsPath) {
    seqan3::sam_file_input alignmentsIn{mappingsInPath, sam_field_ids{}};

    auto &header = alignmentsIn.header();
    std::vector<size_t> ref_lengths{};
    std::ranges::transform(header.ref_id_info, std::back_inserter(ref_lengths),
                           [](auto const &info) { return std::get<0>(info); });

    const std::deque<std::string> &ref_ids = header.ref_ids();

    seqan3::sam_file_output splitsOut{splitsPath, ref_ids, ref_lengths, sam_field_ids{}};
    seqan3::sam_file_output multiSplitsOut{multSplitsPath, ref_ids, ref_lengths, sam_field_ids{}};

    std::string currentRecordId;
    std::vector<SamRecord> currentRecordList;
    std::unordered_set<size_t> invalidRecordHitGroups;
    size_t readsCount = 0;

    for (auto &record : alignmentsIn) {
        readsCount++;

        // This line depends on the SAM file being sorted by read ID
        bool isNewRecordID = !currentRecordId.empty() && (currentRecordId != record.id());

        if (isNewRecordID) {
            // Remove current records that belong to invalid hit groups
            currentRecordList.erase(
                std::remove_if(currentRecordList.begin(), currentRecordList.end(),
                               [&](const SamRecord &record) {
                                   return invalidRecordHitGroups.contains(
                                       record.tags().get<"HI"_tag>());
                               }),
                currentRecordList.end());

            processReadRecords(currentRecordList, splitsOut, multiSplitsOut);
            currentRecordList.clear();
            invalidRecordHitGroups.clear();
        }

        currentRecordId = record.id();

        if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped) ||
            record.sequence().size() < minReadLength || !record.tags().contains("XJ"_tag) ||
            record.tags().get<"XJ"_tag>() < 2 ||
            invalidRecordHitGroups.contains(record.tags().get<"HI"_tag>())) {
            invalidRecordHitGroups.insert(record.tags().get<"HI"_tag>());
            continue;
        }
        currentRecordList.push_back(std::move(record));
    }

    processReadRecords(currentRecordList, splitsOut, multiSplitsOut);

    Logger::log(LogLevel::INFO, "Processed ", readsCount, " reads");
}

void Detect::processReadRecords(const std::vector<SamRecord> &readRecords, auto &splitsOut,
                                auto &multiSplitsOut) {
    if (readRecords.empty()) {
        return;
    }

    const auto splitRecords = getSplitRecords(readRecords);

    if (!splitRecords.has_value()) {
        return;
    }

    // TODO Implement multi split writing
    writeSamFile(splitsOut, splitRecords.value().splitRecords);
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

    auto const addOtherCigar = [&](const auto &cigar) {
        const auto cigarValue = get<0>(cigar);
        endPosRead += cigarValue;
        endPosSplit += cigarValue;
        currentCigar.push_back(cigar);
    };

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

    auto const addDeletionCigar = [&](const auto &cigar) { currentCigar.push_back(cigar); };

    auto const addSoftClipCigar = [&](const auto &cigar) {
        if (!excludeSoftClipping) {
            addOtherCigar(cigar);
        }

        /* If current cigar is empty, we are at the beginning of the read in case of soft
        clipping at the end of the read it is just ignored */
        if (currentCigar.empty()) {
            const auto cigarValue = get<0>(cigar);
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
        referencePosition += cigarValue + endPosRead;
        startPosSplit = endPosSplit;
        startPosRead = endPosRead;
        currentCigar.clear();
    };

    for (const auto &cigar : readRecord.cigar_sequence()) {
        if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
            cigar == 'X'_cigar_operation || cigar == 'I'_cigar_operation) {
            addOtherCigar(cigar);
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

std::optional<EvaluatedSplitRecords> Detect::getSplitRecords(
    const std::vector<SamRecord> &readRecords) {
    std::unordered_map<size_t, std::vector<SamRecord>> recordHitGroups{};
    for (const auto &record : readRecords) {
        recordHitGroups[record.tags().get<"HI"_tag>()].push_back(record);
    }

    std::optional<EvaluatedSplitRecords> bestSplitRecords{};

    const auto insertBestSplitRecords = [&](SplitRecords &splitRecords) {
        const auto evaluatedSplitRecords = EvaluatedSplitRecords::calculateEvaluatedSplitRecords(
            splitRecords, minComplementarity, minFraction, params["nrgmax"].as<double>());

        if (!evaluatedSplitRecords.has_value()) {
            return;
        }

        if (!bestSplitRecords.has_value() ||
            evaluatedSplitRecords.value() > bestSplitRecords.value()) {
            bestSplitRecords.emplace(evaluatedSplitRecords.value());
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
    fs::path tmpInDirPath = mappingsFilePath.parent_path() / "tmpIn";
    fs::path tmpSplitsOutDirPath = splitsFilePath.parent_path() / "tmpSplit";
    fs::path tmpMultiSplitsOutDirPath = splitsFilePath.parent_path() / "tmpMsplit";

    helper::createTmpDir(tmpInDirPath);
    helper::createTmpDir(tmpSplitsOutDirPath);
    helper::createTmpDir(tmpMultiSplitsOutDirPath);

    std::string tmpInPath = (tmpInDirPath / "tmp").string();

    const std::vector<fs::path> splitMappingsFilePaths =
        splitMappingsFile(mappingsFilePath, tmpInPath, mappingRecordsCount);

    std::vector<Detect::ChunkedInOutFilePaths> subFiles;
    for (auto &file : splitMappingsFilePaths) {
        fs::path tmpSplitsOutPath = tmpSplitsOutDirPath / file.filename();
        fs::path tmpMultiSplitsOutPath = tmpMultiSplitsOutDirPath / file.filename();
        subFiles.push_back({file, tmpSplitsOutPath, tmpMultiSplitsOutPath});
    }

    return subFiles;
}

std::vector<fs::path> Detect::splitMappingsFile(const fs::path &mappingsFilePath,
                                                const fs::path &tmpInPath,
                                                const int entries) const {
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

void Detect::createStatisticsFile(const fs::path &splitsFilePath,
                                  const fs::path &multSplitsFilePath,
                                  const fs::path &statsFilePath) const {
    std::ofstream statsFileStream(statsFilePath, std::ios::app);

    if (!statsFileStream.is_open()) {
        throw std::runtime_error("Could not open the stats file.");
    }

    statsFileStream << splitsFilePath.stem().string() << "\t"
                    << helper::countUniqueSamEntries(splitsFilePath) << "\t"
                    << helper::countUniqueSamEntries(multSplitsFilePath) << "\n";
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

    std::vector<ChunkedInOutFilePaths> subFiles =
        prepareInputOutputFiles(mappingRecordsInPath, mergedSplitsOutPath, entries);

#pragma omp parallel for num_threads(threads)
    for (const auto &subFileChunk : subFiles) {
        iterateSortedMappingsFile(subFileChunk.mappingsInPath.string(),
                                  subFileChunk.splitsOutPath.string(),
                                  subFileChunk.multSplitsOutPath.string());
    }

    // Merge results from all processing chunks
    std::vector<fs::path> splitsChunkPaths;
    std::vector<fs::path> multiSplitsChunkPaths;

    for (const auto &subFileChunk : subFiles) {
        splitsChunkPaths.push_back(subFileChunk.splitsOutPath);
        multiSplitsChunkPaths.push_back(subFileChunk.multSplitsOutPath);
    }

    helper::mergeSamFiles(splitsChunkPaths, mergedSplitsOutPath);
    helper::mergeSamFiles(multiSplitsChunkPaths, mergedMultiSplitsOutPath);

    // Remove temporary files
    fs::remove_all(subFiles.front().mappingsInPath.parent_path());
    fs::remove_all(subFiles.front().splitsOutPath.parent_path());
    fs::remove_all(subFiles.front().multSplitsOutPath.parent_path());

    if (params["stats"].as<bool>()) {
        fs::path statsFilePath = output.get<std::string>("stats");
        createStatisticsFile(mergedSplitsOutPath, mergedMultiSplitsOutPath, statsFilePath);
    }
}
