#include "SeqRickshaw.hpp"

SeqRickshaw::SeqRickshaw(const po::variables_map &params) : params(params) {
  readtype = params["readtype"].as<std::string>();

  trimPolyG = params["trimpolyg"].as<bool>();

  adpt5f = params["adpt5f"].as<std::string>();
  adpt5r = params["adpt5r"].as<std::string>();
  adpt3f = params["adpt3f"].as<std::string>();
  adpt3r = params["adpt3r"].as<std::string>();
  missmatchRateTrim = params["mtrim"].as<double>();

  minPhread = params["minqual"].as<int>();
  minLen = params["minlen"].as<int>();
  windowTrimSize = params["wtrim"].as<int>();
  threads = params["threads"].as<int>();
  chunkSize = params["chunksize"].as<int>();

  minOverlapMerge = params["minovl"].as<int>();
  missmatchRateMerge = params["mmerge"].as<double>();
}

void SeqRickshaw::start(pt::ptree sample) {
  if (readtype == "SE") {
    processSingleEnd(sample);
  } else if (readtype == "PE") {
    processPairedEnd(sample);
  } else {
    Logger::log(LogLevel::ERROR, "Invalid readtype specified: " + readtype);
    exit(1);
  }
}

/**
 * Processes a single-end sample.
 *
 * @param sample The sample to be processed.
 */
void SeqRickshaw::processSingleEnd(pt::ptree sample) {
  const pt::ptree input = sample.get_child("input");
  const pt::ptree output = sample.get_child("output");

  const std::string recInPath = input.get<std::string>("forward");
  seqan3::sequence_file_input recIn{recInPath};

  std::string recOutPath = output.get<std::string>("forward");
  seqan3::sequence_file_output recOut{recOutPath};

  const std::string sampleName = std::filesystem::path(recInPath).stem().string();
  Logger::log(LogLevel::INFO, "Started pre-processing " + sampleName + " in single-end mode ...");

  const std::vector<SeqRickshaw::Adapter> adapters5 =
      loadAdapters(adpt5f, SeqRickshaw::TrimConfig::Mode::FIVE_PRIME);
  const std::vector<SeqRickshaw::Adapter> adapters3 =
      loadAdapters(adpt3f, SeqRickshaw::TrimConfig::Mode::THREE_PRIME);

  processSingleEndFileInChunks(recInPath, recOutPath, adapters5, adapters3, chunkSize, threads);

  Logger::log(LogLevel::INFO, "Finished pre-processing " + sampleName + " in single-end mode ...");
}

/**
 * Process the paired-end sample.
 *
 * @param sample The sample to be processed.
 */
void SeqRickshaw::processPairedEnd(pt::ptree sample) {
  const pt::ptree input = sample.get_child("input");
  const pt::ptree output = sample.get_child("output");

  const std::string forwardRecInPath = input.get<std::string>("forward");
  const std::string reverseRecInPath = input.get<std::string>("reverse");

  seqan3::sequence_file_input forwardRecIn{forwardRecInPath};
  seqan3::sequence_file_input reverseRecIn{reverseRecInPath};

  std::string snglFwdOutPath = output.get<std::string>("R1only");
  std::string snglRevOutPath = output.get<std::string>("R2only");

  seqan3::sequence_file_output snglFwdOut{snglFwdOutPath};
  seqan3::sequence_file_output snglRevOut{snglRevOutPath};

  std::string mergeRecOutPath = output.get<std::string>("forward");
  seqan3::sequence_file_output mergeRecOut{mergeRecOutPath};

  std::vector<SeqRickshaw::Adapter> adapters5f =
      loadAdapters(adpt5f, SeqRickshaw::TrimConfig::Mode::FIVE_PRIME);
  std::vector<SeqRickshaw::Adapter> adapters3f =
      loadAdapters(adpt3f, SeqRickshaw::TrimConfig::Mode::THREE_PRIME);
  std::vector<SeqRickshaw::Adapter> adapters5r =
      loadAdapters(adpt5r, SeqRickshaw::TrimConfig::Mode::FIVE_PRIME);
  std::vector<SeqRickshaw::Adapter> adapters3r =
      loadAdapters(adpt3r, SeqRickshaw::TrimConfig::Mode::THREE_PRIME);

  const std::string sampleName = std::filesystem::path(forwardRecInPath).stem().string();
  Logger::log(LogLevel::INFO, "Started pre-processing " + sampleName + " in paired-end mode ...");

  processPairedEndFileInChunks(forwardRecInPath, reverseRecInPath, mergeRecOutPath, snglFwdOutPath,
                               snglRevOutPath, adapters5f, adapters3f, adapters5r, adapters3r,
                               chunkSize, threads);

  Logger::log(LogLevel::INFO, "Finished pre-processing " + sampleName + " in paired-end mode ...");
}

/**
 * @brief Loads adapters from a file or a sequence.
 *
 * This function loads adapters from either a file or a sequence and returns them as a vector of
 * SeqRickshaw::Adapter objects.
 *
 * @param filenameOrSequence The filename or sequence from which to load the adapters.
 * @param trimmingMode The trimming mode to be used for the adapters.
 * @return A vector of SeqRickshaw::Adapter objects containing the loaded adapters.
 */
std::vector<SeqRickshaw::Adapter> SeqRickshaw::loadAdapters(
    std::string const &filenameOrSequence, const SeqRickshaw::TrimConfig::Mode trimmingMode) {
  std::vector<SeqRickshaw::Adapter> adapters;

  if (filenameOrSequence.empty()) return adapters;

  if (std::filesystem::exists(filenameOrSequence)) {
    seqan3::sequence_file_input adapterFile{filenameOrSequence};
    for (auto &record : adapterFile) {
      adapters.push_back(SeqRickshaw::Adapter{record.sequence(), missmatchRateTrim, trimmingMode});
    }
  } else {
    seqan3::dna5_vector adapterSequence = filenameOrSequence |
                                          seqan3::views::char_to<seqan3::dna5> |
                                          seqan3::ranges::to<std::vector>();
    adapters.push_back(SeqRickshaw::Adapter{adapterSequence, missmatchRateTrim, trimmingMode});
  }

  for (auto &adapter : adapters) {
    Logger::log(LogLevel::INFO, "Loaded adapter: " + filenameOrSequence + " with " +
                                    std::to_string(adapter.maxMissmatches) +
                                    " allowed missmatches.");
  }

  return adapters;
}

/**
 * Trims the windowed quality of a given record.
 *
 * @tparam record_type The type of the record.
 * @param record The record to trim.
 */
template <typename record_type>
void SeqRickshaw::trimWindowedQuality(record_type &record) {
  std::size_t trimmingEnd = record.sequence().size();

  while ((trimmingEnd - windowTrimSize) >= windowTrimSize) {
    auto windowQual =
        record.base_qualities() | seqan3::views::slice(trimmingEnd - windowTrimSize, trimmingEnd);

    auto windowPhred =
        windowQual | std::views::transform([](auto quality) { return seqan3::to_phred(quality); });

    auto windowPhredSum = std::accumulate(windowPhred.begin(), windowPhred.end(), 0);
    auto windowMeanPhred = windowPhredSum / std::ranges::size(windowPhred);

    if (windowMeanPhred >= minPhread) {
      break;
    }
    trimmingEnd--;
  }

  record.sequence().erase(record.sequence().begin() + trimmingEnd, record.sequence().end());
  record.base_qualities().erase(record.base_qualities().begin() + trimmingEnd,
                                record.base_qualities().end());
}

/**
 * @brief Merges two records if they have a sufficient overlap.
 *
 * This function takes two records and checks if they have a sufficient overlap based on the
 * specified minimum overlap value. If the overlap is sufficient, the records are merged into a new
 * record and returned as an optional. Otherwise, std::nullopt is returned.
 *
 * @tparam record_type The type of the records.
 * @param record1 The first record to be merged.
 * @param record2 The second record to be merged.
 * @param minOverlap The minimum required overlap between the records.
 * @return An optional containing the merged record if the overlap is sufficient, otherwise
 * std::nullopt.
 */
template <typename record_type>
std::optional<record_type> SeqRickshaw::mergeRecordPair(record_type &record1,
                                                        record_type &record2) {
  seqan3::align_cfg::method_global semiGlobalConfig{
      seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
      seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
      seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
      seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};
  seqan3::align_cfg::scoring_scheme scoringSchemeConfig{
      seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
  seqan3::align_cfg::gap_cost_affine gapSchemeConfig{seqan3::align_cfg::open_score{-2},
                                                     seqan3::align_cfg::extension_score{-4}};

  auto outputConfig = seqan3::align_cfg::output_score{} |
                      seqan3::align_cfg::output_begin_position{} |
                      seqan3::align_cfg::output_alignment{};

  auto alignmentConfig = semiGlobalConfig | scoringSchemeConfig | gapSchemeConfig | outputConfig;

  auto &seq1 = record1.sequence();
  auto seq2ReverseComplement = record2.sequence() | std::views::reverse | seqan3::views::complement;

  for (auto const &result :
       seqan3::align_pairwise(std::tie(seq1, seq2ReverseComplement), alignmentConfig)) {
    const int overlap = seq1.size() - result.sequence1_begin_position();
    if (overlap < minOverlapMerge) continue;

    const int minScore = overlap - (overlap * missmatchRateMerge);
    if (result.score() >= minScore) {
      return constructMergedRecord(record1, record2, overlap);
    }
  }

  return std::nullopt;
}

/**
 * Constructs a merged record by combining two input records with a specified overlap.
 *
 * @tparam record_type The type of the input records.
 * @param record1 The first input record.
 * @param record2 The second input record.
 * @param overlap The length of the overlap between the two records.
 * @return The merged record.
 */
template <typename record_type>
record_type SeqRickshaw::constructMergedRecord(const record_type &record1,
                                               const record_type &record2, const int overlap) {
  auto record2ReverseComplementSequence =
      record2.sequence() | seqan3::views::complement | std::views::reverse;
  auto record2ReverseComplementQualities = record2.base_qualities() | std::views::reverse;

  auto nonOverlapStartSequence =
      record1.sequence() | seqan3::views::slice(0, record1.sequence().size() - overlap);
  auto nonOverlapStartQualities =
      record1.base_qualities() | seqan3::views::slice(0, record1.sequence().size() - overlap);
  auto nonOverlapEndSequence =
      record2ReverseComplementSequence |
      seqan3::views::slice(overlap, record2ReverseComplementSequence.size());
  auto nonOverlapEndQualities =
      record2ReverseComplementQualities |
      seqan3::views::slice(overlap, record2ReverseComplementQualities.size());

  auto overlapSequenceRecord1 =
      record1.sequence() |
      seqan3::views::slice(record1.sequence().size() - overlap, record1.sequence().size());
  auto overlapSequenceRecord2 = record2ReverseComplementSequence | seqan3::views::slice(0, overlap);
  auto overlapQualitiesRecord1 =
      record1.base_qualities() |
      seqan3::views::slice(record1.sequence().size() - overlap, record1.sequence().size());
  auto overlapQualitiesRecord2 =
      record2ReverseComplementQualities | seqan3::views::slice(0, overlap);

  seqan3::dna5_vector overlapSequence;
  std::vector<seqan3::phred42> overlapQualities;

  overlapSequence.reserve(overlap);
  overlapQualities.reserve(overlap);

  for (auto const &[base1, base2, quality1, quality2] :
       seqan3::views::zip(overlapSequenceRecord1, overlapSequenceRecord2, overlapQualitiesRecord1,
                          overlapQualitiesRecord2)) {
    if (base1 == base2) {
      overlapSequence.push_back(base1);
      overlapQualities.push_back(std::max(quality1, quality2));
    } else {
      if (quality1 > quality2) {
        overlapSequence.push_back(base1);
        overlapQualities.push_back(quality1);
      } else {
        overlapSequence.push_back(base2);
        overlapQualities.push_back(quality2);
      }
    }
  }

  seqan3::dna5_vector mergedSequence;
  std::vector<seqan3::phred42> mergedQualities;

  mergedSequence.reserve(nonOverlapStartSequence.size() + overlapSequence.size() +
                         nonOverlapEndSequence.size());
  mergedQualities.reserve(nonOverlapStartQualities.size() + overlapQualities.size() +
                          nonOverlapEndQualities.size());

  mergedSequence.insert(mergedSequence.end(),
                        std::make_move_iterator(nonOverlapStartSequence.begin()),
                        std::make_move_iterator(nonOverlapStartSequence.end()));
  mergedSequence.insert(mergedSequence.end(), std::make_move_iterator(overlapSequence.begin()),
                        std::make_move_iterator(overlapSequence.end()));
  mergedSequence.insert(mergedSequence.end(),
                        std::make_move_iterator(nonOverlapEndSequence.begin()),
                        std::make_move_iterator(nonOverlapEndSequence.end()));

  mergedQualities.insert(mergedQualities.end(),
                         std::make_move_iterator(nonOverlapStartQualities.begin()),
                         std::make_move_iterator(nonOverlapStartQualities.end()));
  mergedQualities.insert(mergedQualities.end(), std::make_move_iterator(overlapQualities.begin()),
                         std::make_move_iterator(overlapQualities.end()));
  mergedQualities.insert(mergedQualities.end(),
                         std::make_move_iterator(nonOverlapEndQualities.begin()),
                         std::make_move_iterator(nonOverlapEndQualities.end()));

  record_type mergedRecord{std::move(mergedSequence), record1.id(), std::move(mergedQualities)};

  return mergedRecord;
}

/**
 * Checks if a record passes the quality and length filters.
 *
 * @param record The record to be checked.
 * @return True if the record passes all the filters, false otherwise.
 */
bool SeqRickshaw::passesFilters(auto &record) {
  // Filter for mean quality
  auto phredQual =
      record.base_qualities() | std::views::transform([](auto q) { return q.to_phred(); });
  double sum = std::accumulate(phredQual.begin(), phredQual.end(), 0);
  double meanPhred = sum / std::ranges::size(phredQual);
  bool passesQual = meanPhred >= minPhread;

  // Filter for length
  bool passesLen = std::ranges::distance(record.sequence()) >= minLen;

  return passesQual && passesLen;
}

/**
 * Trims adapter sequences from the given record.
 *
 * @param adapterSequence The adapter sequence to be trimmed.
 * @param record The record from which the adapter sequence will be trimmed.
 * @param trimmingMode The trimming mode to be applied.
 */
void SeqRickshaw::trimAdapter(const SeqRickshaw::Adapter &adapter, auto &record) {
  seqan3::align_cfg::scoring_scheme scoringSchemeConfig{
      seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
  seqan3::align_cfg::gap_cost_affine gapSchemeConfig{seqan3::align_cfg::open_score{-1},
                                                     seqan3::align_cfg::extension_score{-2}};
  auto outputConfig = seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_end_position{} |
                      seqan3::align_cfg::output_begin_position{};

  auto alignment_config = SeqRickshaw::TrimConfig::alignmentConfigFor(adapter.trimmingMode) |
                          scoringSchemeConfig | gapSchemeConfig | outputConfig;

  auto &seq = record.sequence();
  auto &qual = record.base_qualities();

  for (auto const &result :
       seqan3::align_pairwise(std::tie(adapter.sequence, seq), alignment_config)) {
    int minScore = adapter.sequence.size() - adapter.maxMissmatches;
    if (result.score() >= minScore) {
      if (adapter.trimmingMode == SeqRickshaw::TrimConfig::Mode::FIVE_PRIME) {
        seq.erase(seq.begin(), seq.begin() + result.sequence2_end_position());
        qual.erase(qual.begin(), qual.begin() + result.sequence2_end_position());
      } else {
        seq.erase(seq.begin() + result.sequence2_begin_position(), seq.end());
        qual.erase(qual.begin() + result.sequence2_begin_position(), qual.end());
      }
    }
  }
}

/**
 * @brief Trims trailing polyG bases from the 3' end of a sequence record.
 *
 * This function removes any consecutive 'G' bases from the end of the sequence
 * if they meet the quality threshold. The mean quality threshold is set to a rank of 20
 * using the phred42 scale. If the number of consecutive polyG bases is greater than
 * or equal to 5, they are removed from both the sequence and the base qualities.
 *
 * @tparam record_type The type of the sequence record.
 * @param record The sequence record to trim.
 */
template <typename record_type>
void SeqRickshaw::trim3PolyG(record_type &record) {
  auto &seq = record.sequence();
  auto &qual = record.base_qualities();

  std::vector<seqan3::phred42> qualitiesPhread;
  qualitiesPhread.reserve(qual.size());

  seqan3::phred42 qualityThreshold;
  qualityThreshold.assign_rank(30);

  auto seqIt = seq.rbegin();
  auto qualIt = qual.rbegin();

  auto sufficientMeanQuality = [&]() {
    auto qualities = qualitiesPhread |
                     std::views::transform([](auto quality) { return seqan3::to_phred(quality); });

    auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
    return std::ranges::size(qualities) == 0 || sum / std::ranges::size(qualities) >= 20;
  };

  while (seqIt != seq.rend() && (*seqIt == 'G'_dna5 && sufficientMeanQuality())) {
    qualitiesPhread.push_back(*qualIt);
    ++seqIt;
    ++qualIt;
  }

  std::size_t polyGCount = qualitiesPhread.size();

  if (polyGCount >= 5) {
    seq.erase(seq.end() - polyGCount, seq.end());
    qual.erase(qual.end() - polyGCount, qual.end());
  }
}

/**
 * @brief Returns the semi-global alignment configuration for the given mode, assumes adapter as
 * sequence1.
 *
 * @param mode The mode to use for trimming.
 * @return seqan3::align_cfg::method_global The semi-global alignment configuration.
 */
seqan3::align_cfg::method_global SeqRickshaw::TrimConfig::alignmentConfigFor(
    SeqRickshaw::TrimConfig::Mode mode) {
  seqan3::align_cfg::method_global config;

  switch (mode) {
    case SeqRickshaw::TrimConfig::Mode::FIVE_PRIME:
      config.free_end_gaps_sequence1_leading = true;
      config.free_end_gaps_sequence2_leading = true;
      config.free_end_gaps_sequence1_trailing = false;
      config.free_end_gaps_sequence2_trailing = true;
      break;
    case SeqRickshaw::TrimConfig::Mode::THREE_PRIME:
      config.free_end_gaps_sequence1_leading = false;
      config.free_end_gaps_sequence2_leading = true;
      config.free_end_gaps_sequence1_trailing = true;
      config.free_end_gaps_sequence2_trailing = true;
      break;
  }

  return config;
}

void SeqRickshaw::processSingleEndRecordChunk(SeqRickshaw::SingleEndFastqChunk &chunk,
                                              const std::vector<SeqRickshaw::Adapter> &adapters5,
                                              const std::vector<SeqRickshaw::Adapter> &adapters3) {
  auto it = std::remove_if(chunk.records.begin(), chunk.records.end(), [&](auto &record) {
    if (trimPolyG) trim3PolyG(record);

    if (windowTrimSize > 0) trimWindowedQuality(record);

    for (auto const &adapter : adapters5) trimAdapter(adapter, record);

    for (auto const &adapter : adapters3) trimAdapter(adapter, record);

    return !passesFilters(record);
  });

  chunk.records.erase(it, chunk.records.end());
}

// Function to read a FASTQ file in chunks and process each chunk.
void SeqRickshaw::processSingleEndFileInChunks(std::string const &recInPath, std::string recOutPath,
                                               const std::vector<SeqRickshaw::Adapter> &adapters5,
                                               const std::vector<SeqRickshaw::Adapter> &adapters3,
                                               size_t chunkSize, size_t numThreads) {
  seqan3::sequence_file_input recIn{recInPath};
  seqan3::sequence_file_output recOut{recOutPath};

  // Mutex to synchronize access to shared data structures.
  std::mutex mutex;

  size_t chunkCount = 0;
  // Function to read a chunk from the FASTQ file and process it.
  auto processChunkFunc = [&]() {
    while (true) {
      std::chrono::high_resolution_clock::time_point start =
          std::chrono::high_resolution_clock::now();

      // Read a chunk of records from the FASTQ file.
      SingleEndFastqChunk chunk;
      {
        std::lock_guard<std::mutex> lock(mutex);

        chunk.records.reserve(chunkSize);

        size_t count = 0;

        for (auto &&record : recIn) {
          chunk.records.push_back(std::move(record));

          if (++count == chunkSize) break;
        }

        chunkCount++;
      }

      // Process the chunk.
      processSingleEndRecordChunk(chunk, adapters5, adapters3);
      {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto &&record : chunk.records) {
          recOut.push_back(record);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        Logger::log(LogLevel::INFO,
                    "Processed chunk (" + std::to_string(chunkSize) +
                        " reads). Elapsed time: " + std::to_string(duration.count()) + " seconds.");
      }

      // Check if there are no more records in the file.
      if (recIn.begin() == recIn.end()) break;
    }
  };

  // Create and run threads.
  std::vector<std::thread> threads;
  for (size_t i = 0; i < numThreads; ++i) threads.emplace_back(processChunkFunc);

  // Wait for all threads to finish.
  for (auto &thread : threads) thread.join();
}

// Function to read a FASTQ file in chunks and process each chunk.
void SeqRickshaw::processPairedEndFileInChunks(
    std::string const &recFwdInPath, std::string const &recRevInPath,
    std::string const &mergedOutPath, std::string const &snglFwdOutPath,
    std::string const &snglRevOutPath, const std::vector<SeqRickshaw::Adapter> &adapters5f,
    const std::vector<SeqRickshaw::Adapter> &adapters3f,
    const std::vector<SeqRickshaw::Adapter> &adapters5r,
    const std::vector<SeqRickshaw::Adapter> &adapters3r, size_t chunkSize, size_t numThreads) {
  seqan3::sequence_file_input recFwdIn{recFwdInPath};
  seqan3::sequence_file_input recRevIn{recRevInPath};

  seqan3::sequence_file_output mergedOut{mergedOutPath};
  seqan3::sequence_file_output snglFwdOut{snglFwdOutPath};
  seqan3::sequence_file_output snglRevOut{snglRevOutPath};

  // Mutex to synchronize access to shared data structures.
  std::mutex mutex;

  size_t chunkCount = 0;

  // Function to read a chunk from the FASTQ file and process it.
  auto processChunkFunc = [&]() {
    while (true) {
      std::chrono::high_resolution_clock::time_point start =
          std::chrono::high_resolution_clock::now();

      // Read a chunk of records from the FASTQ file.
      PairedEndFastqChunk chunk;
      {
        std::lock_guard<std::mutex> lock(mutex);

        int count = 0;

        for (auto &&[record1, record2] : seqan3::views::zip(recFwdIn, recRevIn)) {
          chunk.recordsFwd.push_back(std::move(record1));
          chunk.recordsRev.push_back(std::move(record2));

          if (++count == chunkSize) break;
        }

        chunkCount++;
      }

      // Process the chunk.
      processPairedEndRecordChunk(chunk, adapters5f, adapters3f, adapters5r, adapters3r);

      {
        std::lock_guard<std::mutex> lock(mutex);

        for (auto &&record : chunk.recordsMergedRes) {
          mergedOut.push_back(record);
        }
        for (auto &&record : chunk.recordsSnglFwdRes) {
          snglFwdOut.push_back(record);
        }
        for (auto &&record : chunk.recordsSnglRevRes) {
          snglRevOut.push_back(record);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        Logger::log(LogLevel::INFO, "Processed chunk (" + std::to_string(chunkSize) +
                                        " read pairs). Elapsed time: " +
                                        std::to_string(duration.count()) + " seconds.");
      }

      // Check if there are no more records in the file.
      if (recFwdIn.begin() == recFwdIn.end() || recRevIn.begin() == recRevIn.end()) break;
    }
  };

  // Create and run threads.
  std::vector<std::thread> threads;
  for (size_t i = 0; i < numThreads; ++i) threads.emplace_back(processChunkFunc);

  // Wait for all threads to finish.
  for (auto &thread : threads) thread.join();
}

void SeqRickshaw::processPairedEndRecordChunk(SeqRickshaw::PairedEndFastqChunk &chunk,
                                              const std::vector<SeqRickshaw::Adapter> &adapters5f,
                                              const std::vector<SeqRickshaw::Adapter> &adapters3f,
                                              const std::vector<SeqRickshaw::Adapter> &adapters5r,
                                              const std::vector<SeqRickshaw::Adapter> &adapters3r) {
  for (auto &&[record1, record2] : seqan3::views::zip(chunk.recordsFwd, chunk.recordsRev)) {
    if (trimPolyG) {
      trim3PolyG(record1);
      trim3PolyG(record2);
    }

    if (windowTrimSize > 0) {
      trimWindowedQuality(record1);
      trimWindowedQuality(record2);
    }

    for (auto const &adapter : adapters5f) trimAdapter(adapter, record1);

    for (auto const &adapter : adapters3f) trimAdapter(adapter, record1);

    for (auto const &adapter : adapters5r) trimAdapter(adapter, record2);

    for (auto const &adapter : adapters3r) trimAdapter(adapter, record2);

    bool filtFwd = passesFilters(record1);
    bool filtRev = passesFilters(record2);

    if (filtFwd && filtRev) {
      auto mergedRecord = mergeRecordPair(record1, record2);
      if (mergedRecord.has_value()) {
        chunk.recordsMergedRes.push_back(std::move(mergedRecord.value()));
      } else {
        chunk.recordsSnglFwdRes.push_back(std::move(record1));
        chunk.recordsSnglRevRes.push_back(std::move(record2));
      }
    } else {
      if (filtFwd) {
        chunk.recordsSnglFwdRes.push_back(std::move(record1));
      }
      if (filtRev) {
        chunk.recordsSnglRevRes.push_back(std::move(record2));
      }
    }
  }
}