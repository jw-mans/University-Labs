#include "matrix.h"

Matrix::Matrix(int r,int c) : rows(r), cols(c), data(r*c) {}

double& Matrix::operator()(int r,int c) { return data[r*cols+c]; }

double Matrix::operator()(int r,int c) const { return data[r*cols+c]; }