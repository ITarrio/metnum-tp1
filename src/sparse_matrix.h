#ifndef METNUM_TP1_SPARSE_MATRIX_H
#define METNUM_TP1_SPARSE_MATRIX_H

#include <map>
#include <vector>
#include <map>
#include "matrix.h"

using std::vector;
using std::pair;
using std::make_pair;

typedef std::map<size_t, double> bucket;

class sparse_matrix {

    public:
        sparse_matrix(size_t rows, size_t cols) : rows(rows), cols(cols), matrix(cols) {
            trans = false;
        }

        size_t getRows() const {
            return (trans) ? cols : rows;
        }

        size_t getCols() const {
            return (trans) ? rows : cols;
        }

        double get(size_t x, size_t y) const {
            if (trans){
                size_t w = x; x = y; y = w;
            }
            const bucket& x_map = matrix[x];
            auto y_map = x_map.find(y);
            return (y_map != x_map.end()) ? y_map->second : 0;
        }

        void set(size_t x, size_t y, double num) {
            if (trans){
                size_t w = x; x = y; y = w;
            }
            bucket& x_map = matrix[x];
            if (num == 0) {
                x_map.erase(y);
            } else {
                x_map[y] = num;
            }

        }

        void transponse() {
            trans = !trans;
        }

        vector<pair<size_t, size_t>> notZeros() {
            vector<pair<size_t, size_t>> nc;
            for (size_t i = 0; i < cols; ++i) {
                bucket& x = matrix[i];
                for(auto y = x.begin(); y != x.end(); ++y) {
                    pair<size_t, size_t> pair = (trans) ? make_pair(y->first, i) : make_pair(i, y->first);
                    nc.push_back(pair);
                }
            }
            return nc;
        };

        row<double> transposedProductWithVector(const row<double> &b) {
            row<double> result(cols, 0);

            for (size_t i = 0; i < cols; ++i) {
                bucket& column = matrix[i];
                double sum = 0;
                for (auto column_row = column.begin(); column_row != column.end(); ++column_row) {
                    sum += column_row->second * b[column_row->first];
                }
                result[i] = sum;
            }

            return result;
        }

        sparse_matrix transposedByNotTransposedProduct() {
            sparse_matrix result(cols, cols);
            for (size_t i = 0; i < cols; ++i) {
                bucket& column = matrix[i];
                for (size_t j = 0; j <= i; ++j) {
                    bucket& another_column = matrix[j];
                    double sum = 0;
                    for (auto column_row = column.begin(); column_row != column.end(); ++column_row) {
                        auto another_column_row = another_column.find(column_row->first);
                        if(another_column_row != another_column.end()) {
                            sum += column_row->second * another_column_row->second;
                        }
                    }
                    result.set(i,j, sum);
                    if (i != j) {
                        result.set(j, i, sum);
                    }
                }
            }
            return result;
        }


        // Esta funcion asume que la matriz es cuadrada y triangular inferior.
        row<double> solveCholeskySystem(row<double> b){
            // Resuelvo Lz = b
            sparse_matrix transposedLower = fullTranspose();
            size_t z_size = b.size();
            row<double> z(z_size, 0);
            for (size_t i = 0; i < cols; ++i) {
                bucket& column = transposedLower.matrix[i];
                double sumOfRowI = 0;
                double c = 0.0;
                for (auto row_column = column.begin(); row_column != column.end() && row_column->first <= i; row_column++) {
                    double y = (row_column->second * z[row_column->first]) - c;
                    double t = sumOfRowI + y;
                    c = (t - sumOfRowI) - y;
                    sumOfRowI = t;
                }
                if (get(i,i) != 0) {
                    z[i] = (b[i] - sumOfRowI) / get(i,i);
                }
            }
            // Resuelvo L'x = z

            size_t x_size = z.size();
            row<double> x(x_size, 0);

            for (size_t i = 0; i < cols; ++i) {
                bucket& column = matrix[cols-1-i];
                double sumOfRowI = 0;
                double c = 0.0;
                for (auto row_column = column.rbegin(); row_column != column.rend(); row_column++) {
                    double y = (row_column->second * x[row_column->first]) - c;
                    double t = sumOfRowI + y;
                    c = (t - sumOfRowI) - y;
                    sumOfRowI = t;
                }
                if (get(cols-1-i,cols-1-i) != 0) {
                    x[cols-1-i] = (z[cols-1-i] - sumOfRowI) / get(cols-1-i,cols-1-i);
                }
            }
            return x;
        }

        const bucket& column(size_t col) const {
            return matrix[col];
        }

        sparse_matrix fullTranspose() {
            sparse_matrix transposedMatrix(cols,rows);
            for (int i = 0; i < cols; ++i) {
                bucket& column = matrix[i];
                for (auto columnIt = column.begin(); columnIt != column.end() ; columnIt++) {
                    transposedMatrix.set(columnIt->first, i, columnIt->second);
                }
            }
            return transposedMatrix;
        }

        void removeColumn(size_t colNumber) {
            matrix.erase(matrix.begin()+colNumber);
            cols -= 1;
        }

        void resizeHeight(size_t actualRows) {
            rows = actualRows;
        }

private:
        bool trans;
        size_t rows;
        size_t cols;
        vector<bucket> matrix;
};


#endif //METNUM_TP1_SPARSE_MATRIX_H
