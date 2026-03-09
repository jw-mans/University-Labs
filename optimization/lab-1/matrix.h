#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

struct Matrix {
    int rows, cols;
    std::vector<double> data;

    Matrix() : rows(0), cols(0), data() {};
    Matrix(int r,int c);

    double& operator()(int r,int c);
    double operator()(int r,int c) const;
};

#endif