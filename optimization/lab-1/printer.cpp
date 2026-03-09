#include "printer.h"
#include <iostream>
#include <iomanip>

void printTableau(
    const Matrix& tableau,
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels,
    int iteration,
    const std::string& phase
) {

    std::cout << "\n--- " << phase << ", Iteration " << iteration << " ---\n";

    std::cout << "     ";
    for (auto& label : col_labels)
        std::cout << std::setw(8) << label;
    std::cout << "\n";

    for (int i = 0; i < tableau.rows; i++) {

        std::cout << std::setw(4) << row_labels[i] << " ";

        for (int j = 0; j < tableau.cols; j++) {

            double v = tableau(i,j);
            if (std::abs(v) < 1e-10) v = 0;

            std::cout << std::setw(8)
                      << std::fixed
                      << std::setprecision(3)
                      << v;
        }

        std::cout << "\n";
    }

    std::cout << "------------------------------------------\n";
}