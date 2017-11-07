#ifndef METNUM_TP1_CHOLESKY_H
#define METNUM_TP1_CHOLESKY_H

#include "matrix.h"
#include "sparse_matrix.h"

/**
 * Factorize the provided matrix if it is a S.P.D. into an lower triangular matrix.
 *
 * @pre the matrix symmetric and positive-definite
 * @tparam T
 * @param factorized lower triangular matrix aka L
 */
template <typename T>
matrix<T> cholesky_factorization(const matrix<T>& mx) {
    assert(Matrix::isSymmetric(mx));
    size_t n = mx.size();

    matrix<T> L(n, row<T>(n, 0));

    L[0][0] = sqrt(mx[0][0]);
    for (size_t j = 1; j < n; ++j) {
        L[j][0] = mx[j][0] / L[0][0];
    }

    for (size_t i = 1; i < n-1; ++i) {
        T sum_col_i = 0;
        T c = 0.0;
        for (size_t k = 0; k < i; ++k) {
            //Super Kahan!
            T y = pow(L[i][k], 2)-c;
            T t = sum_col_i + y;
            c = (t - sum_col_i) - y;
            sum_col_i = t;
        }
        T lii = mx[i][i] - sum_col_i;
        assert(lii > 0);
        L[i][i] = sqrt(lii);

        for (size_t j = i+1; j < n; ++j) {
            T sum_cols_ij = 0;
            c = 0.0;
            for (size_t k=0; k < i; ++k) {
                //Super Kahan!
                T y = (L[j][k] * L[i][k])-c;
                T t = sum_cols_ij + y;
                c = (t - sum_cols_ij) - y;
                sum_cols_ij = t;
            }

            L[j][i] = (mx[j][i] - sum_cols_ij) / L[i][i];
        }
    }
    T sum_Lnk = 0;
    for (size_t k=0; k < n-1; ++k) {
        sum_Lnk += pow(L[n-1][k],2);
    }
    L[n-1][n-1] = sqrt(mx[n-1][n-1] - sum_Lnk);
    return L;
}



sparse_matrix sparse_cholesky_factorization(sparse_matrix& mx) {
    size_t n = mx.getRows();

    sparse_matrix L(n, n);

    for (size_t j = 0; j < mx.getCols(); ++j) {
        const bucket& column = mx.column(j);
        double sumOfColumn = 0;
        for (auto row = column.begin(); row != column.end() && row->first < j; ++row) {
            size_t i = row->first;
            sumOfColumn += pow(L.get(i,j), 2);
            //std::cout << "L ij hay en i=" << i << "j=" << j << " el total de: " << (L.get(i,j)) << std::endl;
        }
        //std::cout << "mx jj hay en j=" << j << " el total de: " << (mx.get(j,j)) << "y en sum of column: " << sumOfColumn << std::endl;
        if (mx.get(j,j) - sumOfColumn >= 0) {
            //std::cout << "asigno en L" << j << j << " " << sqrt(mx.get(j,j) - sumOfColumn) << " que viene de la resta " << mx.get(j,j) - sumOfColumn << std::endl;
            L.set(j,j, sqrt(mx.get(j,j) - sumOfColumn));
        }

        /*for (size_t i = j + 1; i < mx.getCols(); ++i) {
            double sumOfL = 0;
            for (auto row = column.begin(); row != column.end() && row->first < j; ++row) {
                sumOfL += L.get(row->first, j) * L.get(row->first, i);
            }
            if (L.get(j,j) == 0) std::cout << "NO PODES DIVIDIR POR 0 EN CHOLESKY FORRO" << std::endl;
            if (mx.get(j,i) != 0) std::cout << "MX_ij = " << mx.get(j,i) << " i= " << i << " j=" << j << std::endl;
            L.set(j, i, (1/L.get(j,j))*(mx.get(j,i) - sumOfL));
        }*/
        size_t i = j+1;
        double sumOfL = 0;
        for (auto row = column.begin(); row != column.end() && row->first < j; ++row) {
            sumOfL += L.get(row->first, j) * L.get(row->first, i);
        }
        if (L.get(j,j) != 0) {
            //if (j==1) std::cout << "en " << "i=" << i << " j=" << j << "hay: " << L.get(j,j) << ", " << mx.get(j,i) << ", " << sumOfL << std::endl;
            L.set(j, i, (1/L.get(j,j))*(mx.get(j,i) - sumOfL));
        }

        i = j+340;
        sumOfL = 0;
        for (auto row = column.begin(); row != column.end() && row->first < j; ++row) {
            sumOfL += L.get(row->first, j) * L.get(row->first, i);
        }
        if (L.get(j,j) != 0) {
            //std::cout << "NO PODES DIVIDIR POR 0 EN CHOLESKY FORRO en la col i=" << i << " j=" << j << std::endl;
            //std::cout << " L ji hay en i=" << i << " j=" << j << " el total de: " << (L.get(j,j)) << std::endl;
            //std::cout << "MX ji hay en i=" << i << " j=" << j << " el total de: " << (mx.get(j,i)) << std::endl;
            L.set(j, i, (1/L.get(j,j))*(mx.get(j,i) - sumOfL));
        }

        if (j % 10000 == 0) std::cout << "van j=" << j << "de " << mx.getCols() << std::endl;
    }
    L.transponse();
    return L;
}

#endif //METNUM_TP1_CHOLESKY_H
