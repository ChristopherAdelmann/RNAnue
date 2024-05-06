//
// Created by Richard Albin Schaefer on 1/23/24.
//

#ifndef RNANUE_BASE_HPP
#define RNANUE_BASE_HPP

#include <boost/program_options.hpp>

#include "Data.hpp"
#include "Utility.hpp"

namespace po = boost::program_options;

class Base {
   public:
    Base(po::variables_map _params);
    ~Base() = default;
    Data data;

   private:
    po::variables_map params;
    //        Data data;
};

#endif  // RNANUE_BASE_HPP
