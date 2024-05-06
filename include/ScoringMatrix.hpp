#ifndef SCORINGMATRIX_H
#define SCORINGMATRIX_H

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <bitset>
#include <iostream>
#include <map>

typedef struct {
    size_t width;
    size_t height;
    int *scoring;               // scores
    std::bitset<3> *traceback;  //
} ScoringMatrix;

/* traceback matrix
 * 0 indicates diagonal
 * 1 indicates left
 * 2 indicates right
 */

ScoringMatrix createScoringMatrix(const char *a, const char *b, int matchScore, int mismatchScore,
                                  int gapCost);

// int similarity(const char *a, const char *b);
void freeMatrix(ScoringMatrix *mat);
void printMatrix(ScoringMatrix mat);

// int findMaxElement(int score[3]);

#endif  // SCORINGMATRIX_H
