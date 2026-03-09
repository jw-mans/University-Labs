#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "matrix.h"
#include <vector>
#include <string>
#include <map>

struct SimplexResult {
    Matrix tableau;
    std::vector<std::string> row_labels;
    std::vector<std::string> col_labels;
    std::map<std::string,double> solution;
};

SimplexResult simplexPhase2(
    Matrix tableau,
    std::vector<std::string> row_labels,
    std::vector<std::string> col_labels,
    double eps,
    const std::string& phase
);

#endif