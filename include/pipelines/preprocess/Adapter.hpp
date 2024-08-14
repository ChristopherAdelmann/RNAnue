#pragma once

// Standard
#include <string>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

// Classes
#include "PreprocessParameters.hpp"
#include "TrimConfig.hpp"

namespace pipelines {
namespace preprocess {

template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

struct Adapter {
    const seqan3::dna5_vector sequence;
    const double maxMissMatchFraction;
    const TrimConfig::Mode trimmingMode;

    static std::vector<Adapter> loadAdapters(const AdapterInput &adapterInput,
                                             const double maxMissMatchFraction,
                                             const TrimConfig::Mode trimmingMode) {
        std::vector<Adapter> adapters;

        auto addAdapterFromSequence = [&adapters, &trimmingMode,
                                       maxMissMatchFraction](const seqan3::dna5_vector &sequence) {
            adapters.push_back(Adapter{sequence, maxMissMatchFraction, trimmingMode});
        };

        auto addAdapterFromFile = [&adapters, &trimmingMode,
                                   maxMissMatchFraction](const std::string &filePath) {
            seqan3::sequence_file_input adapterFile{filePath};
            for (const auto &record : adapterFile) {
                adapters.push_back(Adapter{record.sequence(), maxMissMatchFraction, trimmingMode});
            }
        };

        std::visit(
            overloaded{[](const std::monostate &) {}, addAdapterFromSequence, addAdapterFromFile},
            adapterInput);

        for (const auto &adapter : adapters) {
            Logger::log(LogLevel::INFO, "Loaded adapter: ", adapter, " with ",
                        std::to_string(adapter.maxMissMatchFraction),
                        " allowed miss-match fraction.");
        }

        return adapters;
    }
};
}  // namespace preprocess
}  // namespace pipelines
