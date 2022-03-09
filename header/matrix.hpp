#ifndef MATRIX_H
#define MATRIX_H

template<typename field>
class gauss_elimination
{
    public:
    gauss_elimination() {}
    void transpose(field* mat, int dim);
    field reduce(field* mat, int dim);
    bool inverse(field* mat, field* inv, int dim);
    bool solve(field* mat, field* y, int dim);
};

#include "matrix.tpp"

#endif