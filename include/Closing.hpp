#pragma once

// Standard
#include <algorithm>
#include <random>
#include <vector>

// Classes
#include "Logger.hpp"

class Closing {
   private:
    static std::vector<std::string> retrieveQuotes();

   public:
    Closing() = delete;
    ~Closing() = delete;

    static void printQuote();
};