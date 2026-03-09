#include "simplex.h"
#include "printer.h"

#include <iostream>
#include <limits>
#include <cmath>
#include <map>

SimplexResult simplexPhase2(
    Matrix tableau,
    std::vector<std::string> row_labels,
    std::vector<std::string> col_labels,
    double eps,
    const std::string& phase
)
{
    int iteration = 0;

    int rows = tableau.rows;
    int cols = tableau.cols;

    while (true)
    {
        iteration++;

        printTableau(tableau, row_labels, col_labels, iteration, phase);

        bool optimal = true;

        for(int j = 0; j < cols - 1; j++)
        {
            if(tableau(rows - 1, j) > eps)
                optimal = false;
        }

        if(optimal)
        {
            std::cout << phase << " optimal reached\n";
            break;
        }

        int pivot_col = 0;
        double best = tableau(rows - 1, 0);

        for(int j = 1; j < cols - 1; j++)
        {
            if(tableau(rows - 1, j) > best)
            {
                best = tableau(rows - 1, j);
                pivot_col = j;
            }
        }

        int pivot_row = -1;
        double best_ratio = std::numeric_limits<double>::infinity();

        for(int i = 0; i < rows - 1; i++)
        {
            double a = tableau(i, pivot_col);

            if(a > eps)
            {
                double ratio = tableau(i, cols - 1) / a;

                if(ratio < best_ratio)
                {
                    best_ratio = ratio;
                    pivot_row = i;
                }
            }
        }

        if(pivot_row == -1)
        {
            std::cout << "Unbounded\n";
            break;
        }

        double pivot = tableau(pivot_row, pivot_col);

        // нормализация ведущей строки
        for(int j = 0; j < cols; j++)
            tableau(pivot_row, j) /= pivot;

        // обнуление столбца
        for(int i = 0; i < rows; i++)
        {
            if(i == pivot_row) continue;

            double factor = tableau(i, pivot_col);

            for(int j = 0; j < cols; j++)
                tableau(i, j) -= factor * tableau(pivot_row, j);
        }

        row_labels[pivot_row] = col_labels[pivot_col];
    }

    // извлечение решения
    std::map<std::string,double> solution;

    for(int i = 0; i < rows - 1; i++)
    {
        std::string name = row_labels[i];

        if(name != "W")
            solution[name] = tableau(i, cols - 1);
    }

    SimplexResult result;
    result.tableau = tableau;
    result.row_labels = row_labels;
    result.col_labels = col_labels;
    result.solution = solution;

    return result;
}