#include "gtest/gtest.h"
#include "../src/sparse_matrix.h"

TEST(SparseMatrixTest, solveCholeskySystem) {
    int size = 5;
    double value = 1;

    sparse_matrix L(5,5);

    for (size_t i = 0; i < 5; i++) {
        for (size_t j = 0; j <= i; ++j) {
            L.set(i,j,1);
        }
    }

    row<double> b = { 1, 2, 4, 6, 105 };

    row<double> result = L.solveCholeskySystem(b);

    ASSERT_EQ(result[4],  99);
    ASSERT_EQ(result[3], -97);
    ASSERT_EQ(result[2],   0);
    ASSERT_EQ(result[1],  -1);
    ASSERT_EQ(result[0],   0);
}