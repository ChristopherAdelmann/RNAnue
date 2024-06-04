#pragma once

// Boost
#include <boost/program_options.hpp>

// Classes
#include "Constants.hpp"
#include "Data.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;
namespace pi = constants::pipelines;

class Base {
   public:
    explicit Base(po::variables_map _params);
    ~Base() = default;

   private:
    po::variables_map params;
    Data data;
};