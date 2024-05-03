//
// Created by Richard Albin Schaefer on 1/22/24.
//

#ifndef RNANUE_CLOSING_HPP
#define RNANUE_CLOSING_HPP

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

#endif  // RNANUE_CLOSING_HPP
