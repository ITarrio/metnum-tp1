#ifndef METNUM_TP1_GAUSS_H
#define METNUM_TP1_GAUSS_H

#include "matrix.h"

template <typename T>
struct PLUMatrix {

    PLUMatrix( matrix<T> P , matrix<T> L, matrix<T> U) : P(P), L(L), U(U) {}
    PLUMatrix(size_t matrixSize) {
        P = Matrix::identityMatrix<T>((int) matrixSize);
        L = Matrix::zeroMatrix<T>((int) matrixSize);
        U = Matrix::zeroMatrix<T>((int) matrixSize);
    };

    matrix<T> P;
    matrix<T> L;
    matrix<T> U;

};

template <typename T>
PLUMatrix<T> pluFactorization(const matrix<T> &mx);

/**
 * Searches for a row with a non-zero value in a particular column,
 * and swaps it into the row of the same number if needed.
 *
 * @param m the matrix to swap rows
 * @param col the column where the non-zero value should be
 * @return whether such a row was found
 */
template <typename T>
size_t find_swap_rows(matrix<T>& m, size_t col);

template <typename T>
PLUMatrix<T> pluFactorization(const matrix<T> &original_matrix) {
    matrix<T> mx = original_matrix;
    PLUMatrix<T> plu((int) mx.size());
    for (size_t i = 0; i < mx.size() - 1; ++i) {
        size_t destinationRow = find_swap_rows(mx, i); //this should be done only if mx[i][i] is zero?
        if (destinationRow != i) {
            matrix<T> identity(Matrix::identityMatrix<T>((int) mx.size()));
            Matrix::swap_rows(identity, i,destinationRow);
            Matrix::swap_rows(plu.L, i,destinationRow);
            plu.P = Matrix::dotProduct(identity, plu.P);
        }
        if (mx[i][i] != 0) {
            for (size_t j = i+1; j < mx.size(); ++j) {
                plu.L[j][i] = (mx[j][i] / mx[i][i]);
                mx[j][i] = 0;
                for (size_t k = 0; k < mx.size(); k++) {
                    if (k != i) {
                        mx[j][k] -= plu.L[j][i] * mx[i][k];
                    }
                }
            }
        }
    }
    plu.L = Matrix::sum(plu.L, Matrix::identityMatrix<T>((int) plu.L.size()));
    plu.U = mx;
    return plu;
}

template<typename T>
size_t find_swap_rows(matrix<T>& m, size_t col) {
    size_t max = col;
    for (size_t i = col + 1; i < m.size(); ++i) {
        if((m[i][col] > m[max][col] && m[i][col] != 0) || (m[max][col] == 0 && m[i][col] != 0)) {
            max = i;
        }
    }
    if(max != col) {
        std::swap(m[col], m[max]);
    }
    return max;
}

#endif //METNUM_TP1_GAUSS_H
