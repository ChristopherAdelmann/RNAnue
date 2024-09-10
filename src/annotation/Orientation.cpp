#include "Orientation.hpp"

namespace annotation {

std::istream &operator>>(std::istream &in, Orientation &orientation) {
    std::string token;
    in >> token;
    if (token == "both") {
        orientation = Orientation::BOTH;
    } else if (token == "opposite") {
        orientation = Orientation::OPPOSITE;
    } else if (token == "same") {
        orientation = Orientation::SAME;
    } else {
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

}  // namespace annotation
