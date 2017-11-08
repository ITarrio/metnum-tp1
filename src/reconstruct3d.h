#ifndef METNUM_TP1_RECONSTRUCT3D_H
#define METNUM_TP1_RECONSTRUCT3D_H


#include "matrix.h"
#include "plu.h"
#include "cholesky.h"
#include "calibration.h"
#include "sparse_matrix.h"

struct options {
    options(bool maskOptimization) : maskOptimization(maskOptimization) {}
    bool maskOptimization;
};

matrix<double> sourceOfLightMatrix(const direction &s1, const direction &s2, const direction &s3) {
    return {
            {s1.x, s1.y, s1.z},
            {s2.x, s2.y, s2.z},
            {s3.x, s3.y, s3.z}
    };
}
/**
 * (5):
 * for each pixel
 *  /s1x s2x s3x \  /mx \    /I1 \
 * | s1y s2y s3y | | my | = | I2 |
 * \ s1z s2z s3z/  \ mz/    \ I3/
 *
 *
 * (6)
 * agarrar el vector de m y a cada componente la dividis por la norma m
 *
 * Construction of the normal field, for each pixel
 *
 * @param i1, i2, i3 the matrix with the brightness on each pixel for the 3 images
 * @param s1, s2, s3 the source of light for the images
 * @return the normal field of the image
 */
matrix<row<double>> normalField(const matrix<double> &i1, const matrix<double> &i2, const matrix<double> &i3,
                                const direction &s1, const direction &s2, const direction &s3, options opts) {
    size_t height = i1.size(), width = i1[0].size();
    PLUMatrix<double> sm = pluFactorization(sourceOfLightMatrix(s1, s2, s3));
    matrix<row<double>> normal;

    for (size_t i = 0; i < height; i++) {
        row<row<double>> r;
        for (size_t j = 0; j < width; j++) {
            if (opts.maskOptimization || (i1[i][j] == 0 && i2[i][j] == 0 && i3[i][j] == 0)) {
                r.push_back({0, 0, 0}); //TODO: revisen esto
            } else {
                //(5)
                row<double> m = Matrix::solvePLUSystem(sm.P, sm.L, sm.U, {i1[i][j], i2[i][j], i3[i][j]});
                //(6)
                double mNorm = Matrix::twoNorm(m);
                //if (m[0] != 0 && m[1] != 0 && m[2] != 0 && mNorm != 0) std::cout << "no es 0 i=" << i << " j=" << j << std::endl;
                if(mNorm != 0) {
                    m[0] /= mNorm;
                    m[1] /= mNorm;
                    m[2] /= mNorm;
                }
                r.push_back(m);
            }
        }
        normal.push_back(r);
    }

    return normal;
}

sparse_matrix calculateM(const matrix<row<double>> &n) {
    size_t height = n.size();
    size_t width = n[0].size();
    size_t n_size = width * height;
    sparse_matrix M(2*n_size, n_size);
    size_t n_i = 0;

    // En la construccion de la M hay que salvar los bordes que no tienen posicion borde+1.
    for (size_t x = 0; x < width; x++) {
        for (size_t y = 0; y < height - 1; y++) {
            M.set(n_i, n_i, -n[y][x][2]); //le pongo -nz
            M.set(n_i + 1, n_i, n[y][x][2]); //le pongo nz
            n_i++;
        }
        n_i++;
    }


    n_i = 0; // arrancamos de nuevo de la columa 0, pero usamos un offset de filas

    for (size_t x = 0; x < width - 1; x++) {
        for (size_t y = 0; y < height; y++) {
            M.set(n_i, n_i + n_size, -n[y][x][2]); //le pongo -nz
            //if(n_i + height < n_size) {
            M.set(n_i + height, n_i + n_size, n[y][x][2]); //le pongo nz
            //}
            n_i++;
        }
    }
    return M;
}

/*
sparse_matrix calculateM(const matrix<row<double>> &n) {
    size_t height = n.size(), width = n[0].size(), n_size = width * height;
    sparse_matrix M(2*n_size, n_size);
    for (size_t y = 0; y < height-1; y++) {
        for (size_t x = 0; x < width - 1; x++) {
            M.set(x, y, -n[y][x][2]);
            M.set(x+1, y, n[y][x][2]);
        }
    }
    for (size_t y = 0; y < height-1; y++) {
        for (size_t x = 0; x < width - 1; x++) {
            M.set(x, y, -n[y][x][2]);
            M.set(x+height, y, n[y][x][2]);
        }
    }
    return M;
}*/

vector<double> calculateV(const matrix<row<double>> &n) {
    size_t height = n.size(), width = n[0].size();
    vector<double> v;
    for (size_t x = 0; x < width; x++) {
        for (size_t y = 0; y < height; y++) {
            v.push_back(-n[y][x][1]); //le pongo -ny
        }
    }
    for (size_t x = 0; x < width; x++) {
        for (size_t y = 0; y < height; y++) {
            v.push_back(-n[y][x][0]); //le pongo -nx
        }
    }
    return v;
}

matrix<double> solutionToMatrix(row<double> &z, size_t height, size_t width) {
    matrix<double> result(height, row<double>(width));
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            result[i][j] = z[i + j * height];
        }
    }
    return result;
}

//Aqui viene lo bueno jovenes, I cho cho choleskyou
matrix<double> findDepth(const matrix<row<double>> &normalField) {
    std::cout << "Antes de calcular M" << std::endl;
    sparse_matrix m = calculateM(normalField);

    //std::cout << "Despues de calcular M y antes de V" << std::endl;
    vector<double> v = calculateV(normalField);

    std::cout << "Despues de calcular V y antes de calcular MtM" << std::endl;
    // Calcular Mtraspuesta*M : Step 14a
    sparse_matrix mtm = m.transposedByNotTransposedProduct();

    std::cout << "Despues de calcular MtM y antes de calcular Mtv" << std::endl;
    // Calcular Mtraspuesta*V : Step 14b
    row<double> b = m.transposedProductWithVector(v);

    std::cout << "Despues de calcular Mtv" << std::endl;
    // Choleskiar Az=b : Step 15

    std::cout << "Antes de factorizar" << std::endl;
    sparse_matrix L = sparse_cholesky_factorization(mtm);

    std::cout << "Despues de factorizar y antes de resolver" << std::endl;
    row<double> z = L.solveCholeskySystem(b);

    std::cout << "Despues de choleskiar y antes de transformar vector a matriz" << std::endl;
    // z tiene forma (z11,z12,...,z1h,z21,z22,...,z2h,z31,...,zwh)
    size_t height = normalField.size();
    size_t width = normalField[0].size();
    matrix<double> d = solutionToMatrix(z, height, width);

    std::cout << "Despues de transformar vector a matriz" << std::endl;
    return d;
}

#endif //METNUM_TP1_RECONSTRUCT3D_H
