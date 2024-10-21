#pragma once

// Standard
#include <istream>

namespace annotation {

enum Orientation { SAME, OPPOSITE, BOTH };
auto operator>>(std::istream& input, Orientation& orientation) -> std::istream&;

}  // namespace annotation
