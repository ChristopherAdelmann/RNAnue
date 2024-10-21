#include "Orientation.hpp"

namespace annotation {

auto operator>>(std::istream &input, Orientation &orientation) -> std::istream & {
    std::string token;
    input >> token;

    if (token == "both") {
        orientation = Orientation::BOTH;
    } else if (token == "opposite") {
        orientation = Orientation::OPPOSITE;
    } else if (token == "same") {
        orientation = Orientation::SAME;
    } else {
        input.setstate(std::ios_base::failbit);
    }
    return input;
}

}  // namespace annotation
