#pragma once

// Standard
#include <string>
#include <vector>

class Closing {
   private:
    static auto retrieveQuotes() -> std::vector<std::string>;

   public:
    static void printQuote();
};
