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
                x_map.insert(make_pair(y, num));
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
            size_t elementos_totales = 0;
            for (size_t i = 0; i < cols; ++i) {
                size_t elementos = 0;
                bucket& column = matrix[i];
                /*for (size_t j = 0; j <= i; ++j) {
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
                    if (sum != 0) {
                        elementos++;
                        elementos_totales++;
                    }
                }
                if (elementos > 3) std::cout << "total elementos > 3 en col i: " << elementos << ", " << i << std::endl;
                */
                if (i == 0) {
                    size_t j = 0;
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
                    if (sum != 0) {
                        elementos++;
                        elementos_totales++;
                    }
                } else if (i < 340) {
                    size_t j1 = 0;
                    size_t j2 = i-1;
                    bucket& another_column1 = matrix[j1];
                    bucket& another_column2 = matrix[j2];
                    double sum = 0;
                    for (auto column_row = column.begin(); column_row != column.end(); ++column_row) {
                        auto another_column_row = another_column1.find(column_row->first);
                        if(another_column_row != another_column1.end()) {
                            sum += column_row->second * another_column_row->second;
                        }
                    }
                    result.set(i,j1, sum);
                    if (i != j1) {
                        result.set(j1, i, sum);
                    }
                    if (sum != 0) {
                        elementos++;
                        elementos_totales++;
                    }
                    sum = 0;
                    for (auto column_row = column.begin(); column_row != column.end(); ++column_row) {
                        auto another_column_row = another_column2.find(column_row->first);
                        if(another_column_row != another_column2.end()) {
                            sum += column_row->second * another_column_row->second;
                        }
                    }
                    result.set(i,j2, sum);
                    if (i != j2) {
                        result.set(j2, i, sum);
                    }
                    if (sum != 0) {
                        elementos++;
                        elementos_totales++;
                    }
                } else {
                    for (size_t j = i-340; j <= i; ++j) {
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
                        if (sum != 0) {
                            elementos++;
                            elementos_totales++;
                        }
                    }
                }
                if (i % 10000 == 0) std::cout << i << std::endl;
            }
            if (elementos_totales > 0) std::cout << "total elementos en MtM" << elementos_totales << std::endl;
            return result;
        }


        // Esta funcion asume que la matriz es cuadrada y triangular inferior.
        row<double> solveCholeskySystem(row<double> b){
            // Resuelvo Lz = b
            size_t z_size = b.size();
            row<double> z(z_size, 0);
            for (size_t i = 0; i < cols; ++i) {
                bucket& column = matrix[i];
                double sumOfRowI = 0;
                double c = 0.0;
                for (auto row_column = column.begin(); row_column != column.end() && row_column->first <= i; row_column++) {
                    double y = (get(i,row_column->first) * z[row_column->first]) - c;
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

            trans = !trans;
            for (size_t i = cols-1; i > 0; --i) {
                bucket& column = matrix[i];
                double sumOfRowI = 0;
                double c = 0.0;
                for (size_t j = cols-1; j >= i; j--) {
                    double y = (get(i,j) * x[j]) - c;
                    double t = sumOfRowI + y;
                    c = (t - sumOfRowI) - y;
                    sumOfRowI = t;
                }
                if (get(i,i) != 0) {
                    x[i] = (z[i] - sumOfRowI) / get(i,i);
                }
            }
            trans = !trans;
            return x;
        }

        const bucket& column(size_t col) const {
            return matrix[col];
        }

private:
        bool trans;
        size_t rows;
        size_t cols;
        vector<bucket> matrix;
};


#endif //METNUM_TP1_SPARSE_MATRIX_H
