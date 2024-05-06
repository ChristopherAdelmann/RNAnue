#ifndef TRACEBACK_H
#define TRACEBACK_H

#include <vector>

#include "ScoringMatrix.hpp"

typedef struct {
    std::string a;
    std::string b;
    int length;
    int matches;
    int score;
    double cmpl;   // complementarity
    double ratio;  // sitelenratio
} TracebackResult;

std::vector<TracebackResult> traceback(ScoringMatrix matrix, const char *a, const char *b);
void free_traceback(TracebackResult trc);

#endif  // TRACEBACK_H
