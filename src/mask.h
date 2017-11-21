#ifndef METNUM_TP1_MASK_H
#define METNUM_TP1_MASK_H

#include <iostream>
#include "matrix.h"
#include "utils.h"

using namespace std;

namespace Mask {

    struct rect {
        size_t top;
        size_t left;
        size_t bottom;
        size_t right;
    };

    struct mask {
        matrix<bool> img;
        rect clip;
    };

    mask load_mask(const std::string& fname);

    template<typename T>
    matrix<T> apply_clip(const matrix<T>& mx, const mask& msk) {
        const rect& clip = msk.clip;
        size_t height = 1 + clip.bottom - clip.top;
        size_t width = 1 + clip.right - clip.left;
        matrix<T> res;
        for (size_t i = 0; i < height; ++i) {
            const row<T>& r = mx[i + clip.top];
            auto start = r.begin() + clip.left;
            auto end = start + width;
            res.push_back(row<T>(start, end));
        }
        return res;
    }

    template<typename T>
    matrix<T> restore_clip(const matrix<T>& mx, const mask& msk, const T& def = 0) {
        const rect& clip = msk.clip;
        size_t height = msk.img.size();
        size_t width = msk.img[0].size();
        matrix<T> res(clip.top, row<T>(width, def));
        for (const row<T>& r : mx) {
            row<T> new_r(clip.left);
            new_r.insert(new_r.end(), r.begin(), r.end());
            new_r.resize(width, def);
            res.push_back(new_r);
        }
        res.resize(height, row<T>(width, def));
        return res;
    }

    /**
     * Apply a mask to an image, returning a matrix adjusted to the mask size. The
     * elements outside the mask are filled with zeros.
     *
     * @pre The image and mask must have same dimensions.
     * @tparam matrix<T>
     * @param img the image matrix, mask the mask matrix
     */
    template<typename T>
    matrix<T> apply_mask(const matrix<T> &img, const mask& msk) {
        assert(Matrix::rows(img) == Matrix::rows(msk.img)
               && Matrix::columns(img) == Matrix::columns(msk.img));

        matrix<T> result(img);
        for (size_t i = 0; i < result.size(); i++) {
            for (size_t j = 0; j < result[0].size(); j++) {
                if (!msk.img[i][j]) {
                    result[i][j] = 0;
                }
            }
        }

        return result;
    }
}
#endif //METNUM_TP1_MASK_H
