#pragma once

// Standard
#include <istream>

namespace annotation {

enum Orientation { SAME, OPPOSITE, BOTH };
std::istream& operator>>(std::istream& in, Orientation& orientation);

}  // namespace annotation
