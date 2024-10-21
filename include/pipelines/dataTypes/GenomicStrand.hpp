#pragma once

namespace dataTypes {
enum Strand : char { FORWARD = '+', REVERSE = '-' };

inline auto operator!(Strand strand) -> dataTypes::Strand {
    return strand == dataTypes::FORWARD ? dataTypes::REVERSE : dataTypes::FORWARD;
};

}  // namespace dataTypes
