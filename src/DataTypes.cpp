#include "DataTypes.hpp"

namespace dtp {

std::optional<int32_t> recordEndPosition(const SamRecord &record) {
    const auto start = record.reference_position();

    if (!start.has_value()) {
        return std::nullopt;
    }

    int32_t end = start.value();

    for (const auto &cigar : record.cigar_sequence()) {
        if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
            cigar == 'D'_cigar_operation || cigar == 'N'_cigar_operation ||
            cigar == 'X'_cigar_operation) {
            end += get<0>(cigar) - 1;
        }
    }
    return end;
}

GenomicRegion::GenomicRegion(const std::string &referenceID, int32_t startPosition,
                             int32_t endPosition, std::optional<Strand> strand)
    : referenceID(referenceID),
      startPosition(startPosition),
      endPosition(endPosition),
      strand(strand) {}

std::optional<GenomicRegion> GenomicRegion::fromSamRecord(
    const SamRecord &record, const std::deque<std::string> &referenceIDs) {
    const auto start = record.reference_position();
    const auto end = recordEndPosition(record);

    if (!start.has_value() || !end.has_value()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const Strand strand{isReverseStrand ? Strand::REVERSE : Strand::FORWARD};

    return GenomicRegion{referenceIDs[record.reference_id().value()], start.value(), end.value(),
                         strand};
}

}  // namespace dtp