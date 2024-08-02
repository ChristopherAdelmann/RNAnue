#pragma once

// Standard
#include <istream>

namespace Annotation {
enum Orientation { SAME, OPPOSITE, BOTH };
std::istream& operator>>(std::istream& in, Orientation& orientation);
}  // namespace Annotation
