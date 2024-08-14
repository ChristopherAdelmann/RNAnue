#include "Preprocess.hpp"

#include <variant>

#include "Logger.hpp"
#include "PairedRecordMerger.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"
#include "boost/filesystem/path.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/io/sequence_file/input.hpp"

namespace pipelines {
namespace preprocess {
Preprocess::Preprocess(const PreprocessParameters &params) : parameters(params) {}

void Preprocess::start(const PreprocessData &data) {
    Logger::log(LogLevel::INFO, "Processing treatment data");

    for (const auto &sample : data.treatmentSamples) {
        processSample(sample);
    }

    if (!data.controlSamples.has_value()) {
        Logger::log(LogLevel::INFO, "No control samples provided");
        return;
    }

    Logger::log(LogLevel::INFO, "Processing control data");

    for (const auto &sample : data.controlSamples.value()) {
        processSample(sample);
    }
}

void Preprocess::processSample(const PreprocessSampleType &sample) {
    std::visit(
        overloaded{[this](const PreprocessSampleSingle &sample) { processSingleEnd(sample); },
                   [this](const PreprocessSamplePaired &sample) { processPairedEnd(sample); }},
        sample);
}

/**
 * Processes a single-end sample.
 *
 * @param sample The sample to be processed.
 */
void Preprocess::processSingleEnd(const PreprocessSampleSingle &sample) {
    Logger::log(LogLevel::INFO, "Processing sample (single-end): ", sample.input.sampleName);

    const std::vector<Adapter> adapters5 =
        Adapter::loadAdapters(parameters.adapter5Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::FIVE_PRIME);
    const std::vector<Adapter> adapters3 =
        Adapter::loadAdapters(parameters.adapter3Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::THREE_PRIME);

    processSingleEndFileInChunks(sample.input.inputFastqPath, sample.output.outputFastqPath,
                                 adapters5, adapters3, chunkSize, threadCount);

    Logger::log(LogLevel::INFO, "Finished processing sample: ", sample.input.sampleName);
}

/**
 * Process a paired-end sample.
 *
 * @param sample The sample to be processed.
 */
void Preprocess::processPairedEnd(const PreprocessSamplePaired &sample) {
    Logger::log(LogLevel::INFO, "Processing sample (paired-end): ", sample.input.sampleName);

    std::vector<Adapter> adapters5f =
        Adapter::loadAdapters(parameters.adapter5Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::FIVE_PRIME);
    std::vector<Adapter> adapters3f =
        Adapter::loadAdapters(parameters.adapter3Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::THREE_PRIME);
    std::vector<Adapter> adapters5r =
        Adapter::loadAdapters(parameters.adapter5Reverse, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::FIVE_PRIME);
    std::vector<Adapter> adapters3r =
        Adapter::loadAdapters(parameters.adapter3Reverse, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::THREE_PRIME);

    processPairedEndFileInChunks(
        sample.input.inputForwardFastqPath, sample.input.inputReverseFastqPath,
        sample.output.outputMergedFastqPath, sample.output.outputSingletonForwardFastqPath,
        sample.output.outputSingletonReverseFastqPath, adapters5f, adapters3f, adapters5r,
        adapters3r, chunkSize, threadCount);

    Logger::log(LogLevel::INFO, "Finished processing sample: ", sample.input.sampleName);
}

/**
 * Checks if a record passes the quality and length filters.
 *
 * @param record The record to be checked.
 * @return True if the record passes all the filters, false otherwise.
 */
bool Preprocess::passesFilters(const auto &record) {
    // Filter for mean quality
    const auto phredQual =
        record.base_qualities() | std::views::transform([](auto q) { return q.to_phred(); });
    const double sum = std::accumulate(phredQual.begin(), phredQual.end(), 0);
    const double meanPhred = sum / std::ranges::size(phredQual);
    const bool passesQual = meanPhred >= minMeanPhread;

    // Filter for length
    const bool passesLen = record.sequence().size() >= minLen;

    return passesQual && passesLen;
}

void Preprocess::processSingleEndRecordChunk(Preprocess::SingleEndFastqChunk &chunk,
                                             const std::vector<Adapter> &adapters5,
                                             const std::vector<Adapter> &adapters3) {
    auto it = std::remove_if(chunk.records.begin(), chunk.records.end(), [&](auto &record) {
        if (trimPolyG) RecordTrimmer::trim3PolyG(record);

        if (windowTrimSize > 0)
            RecordTrimmer::trimWindowedQuality(record, parameters.windowTrimmingSize,
                                               parameters.minMeanWindowQuality);

        for (auto const &adapter : adapters5)
            RecordTrimmer::trimAdapter(adapter, record, parameters.minOverlapTrimming);

        for (auto const &adapter : adapters3)
            RecordTrimmer::trimAdapter(adapter, record, parameters.minOverlapTrimming);

        return !passesFilters(record);
    });

    chunk.records.erase(it, chunk.records.end());
}

// Function to read a FASTQ file in chunks and process each chunk.
void Preprocess::processSingleEndFileInChunks(std::string const &recInPath, std::string recOutPath,
                                              const std::vector<Adapter> &adapters5,
                                              const std::vector<Adapter> &adapters3,
                                              size_t chunkSize, size_t numThreads) {
    seqan3::sequence_file_input recIn{recInPath};
    seqan3::sequence_file_output recOut{recOutPath};

    // Mutex to synchronize access to shared data structures.
    std::mutex mutex;

    size_t chunkCount = 0;
    // Function to read a chunk from the FASTQ file and process it.
    auto processChunkFunc = [&]() {
        while (true) {
            const std::chrono::high_resolution_clock::time_point start =
                std::chrono::high_resolution_clock::now();

            // Read a chunk of records from the FASTQ file.
            SingleEndFastqChunk chunk;
            {
                std::lock_guard<std::mutex> lock(mutex);

                chunk.records.reserve(chunkSize);

                size_t count = 0;

                for (auto &&record : recIn) {
                    chunk.records.push_back(record);

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

                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> duration = end - start;

                Logger::log(LogLevel::INFO, "Processed chunk (" + std::to_string(chunkSize) +
                                                " reads). Elapsed time: " +
                                                std::to_string(duration.count()) + " seconds.");
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
void Preprocess::processPairedEndFileInChunks(
    const std::string &recFwdInPath, std::string const &recRevInPath,
    const std::string &mergedOutPath, std::string const &snglFwdOutPath,
    const std::string &snglRevOutPath, const std::vector<Adapter> &adapters5f,
    const std::vector<Adapter> &adapters3f, const std::vector<Adapter> &adapters5r,
    const std::vector<Adapter> &adapters3r, size_t chunkSize, size_t numThreads) {
    seqan3::sequence_file_input recFwdIn{recFwdInPath};
    seqan3::sequence_file_input recRevIn{recRevInPath};

    seqan3::sequence_file_output mergedOut{mergedOutPath};
    seqan3::sequence_file_output snglFwdOut{snglFwdOutPath};
    seqan3::sequence_file_output snglRevOut{snglRevOutPath};

    // Mutex to synchronize access to shared data structures.
    std::mutex mutex;

    size_t chunkCount = 0;

    // Function to read a chunk from the FASTQ file and process it.
    const auto processChunkFunc = [&]() {
        while (true) {
            const std::chrono::high_resolution_clock::time_point start =
                std::chrono::high_resolution_clock::now();

            // Read a chunk of records from the FASTQ file.
            PairedEndFastqChunk chunk;
            {
                std::lock_guard<std::mutex> lock(mutex);

                chunk.recordsFwd.reserve(chunkSize);
                chunk.recordsRev.reserve(chunkSize);

                std::size_t count = 0;

                for (auto &&[record1, record2] : seqan3::views::zip(recFwdIn, recRevIn)) {
                    chunk.recordsFwd.push_back(record1);
                    chunk.recordsRev.push_back(record2);

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

                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> duration = end - start;

                Logger::log(LogLevel::INFO, "Processed chunk (up to " + std::to_string(chunkSize) +
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

void Preprocess::processPairedEndRecordChunk(Preprocess::PairedEndFastqChunk &chunk,
                                             const std::vector<Adapter> &adapters5f,
                                             const std::vector<Adapter> &adapters3f,
                                             const std::vector<Adapter> &adapters5r,
                                             const std::vector<Adapter> &adapters3r) {
    for (auto &&[record1, record2] : seqan3::views::zip(chunk.recordsFwd, chunk.recordsRev)) {
        if (trimPolyG) {
            RecordTrimmer::trim3PolyG(record1);
            RecordTrimmer::trim3PolyG(record2);
        }

        if (windowTrimSize > 0) {
            RecordTrimmer::trimWindowedQuality(record1, parameters.windowTrimmingSize,
                                               parameters.minMeanWindowQuality);
            RecordTrimmer::trimWindowedQuality(record2, parameters.windowTrimmingSize,
                                               parameters.minMeanWindowQuality);
        }

        for (auto const &adapter : adapters5f)
            RecordTrimmer::trimAdapter(adapter, record1, parameters.minOverlapTrimming);

        for (auto const &adapter : adapters3f)
            RecordTrimmer::trimAdapter(adapter, record1, parameters.minOverlapTrimming);

        for (auto const &adapter : adapters5r)
            RecordTrimmer::trimAdapter(adapter, record2, parameters.minOverlapTrimming);

        for (auto const &adapter : adapters3r)
            RecordTrimmer::trimAdapter(adapter, record2, parameters.minOverlapTrimming);

        const bool filtFwd = passesFilters(record1);
        const bool filtRev = passesFilters(record2);

        if (filtFwd && filtRev) {
            const auto mergedRecord =
                PairedRecordMerger::mergeRecordPair(record1, record2, parameters.minOverlapMerging,
                                                    parameters.maxMissMatchFractionMerging);

            if (!mergedRecord.has_value()) {
                chunk.recordsSnglFwdRes.push_back(std::move(record1));
                chunk.recordsSnglRevRes.push_back(std::move(record2));
                continue;
            }

            const auto &mergedRecordValue = mergedRecord.value();

            if (passesFilters(mergedRecordValue)) {
                chunk.recordsMergedRes.push_back(std::move(mergedRecord.value()));
                continue;
            }

            continue;
        }

        if (filtFwd) {
            chunk.recordsSnglFwdRes.push_back(std::move(record1));
        }
        if (filtRev) {
            chunk.recordsSnglRevRes.push_back(std::move(record2));
        }
    }
}

}  // namespace preprocess
}  // namespace pipelines
