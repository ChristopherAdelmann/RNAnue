#pragma once

// Standard
#include <string>
#include <vector>

class Closing {
   private:
    std::vector<std::string> quotes;

   public:
    Closing();
    ~Closing();

    std::vector<std::string> retrieveQuotes();
    void printQuote(std::ostream& out);
};