#ifndef PRINTER_H
#define PRINTER_H

#include "matrix.h"
#include <vector>
#include <string>

void printTableau(
    const Matrix& tableau,
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels,
    int iteration,
    const std::string& phase
);

#endif