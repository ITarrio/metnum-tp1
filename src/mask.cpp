#include "mask.h"

using namespace std;

namespace Mask {

    mask load_mask(const std::string& fname) {
        matrix<double> gray = Utils::loadGrayImage(fname);
        size_t height = gray.size();
        size_t width = gray[0].size();
        mask result {
                {height, row<bool>(width, true)},
                {height,width,0,0}
        };
        matrix<bool>& map = result.img;
        rect& clip = result.clip;
        for (size_t i = 0; i < height; ++i) {
            for (size_t j = 0; j < width; ++j) {
                map[i][j] = gray[i][j] > 250;
                if(map[i][j]) {
                    if(clip.left > j) {
                        clip.left = j;
                    }
                    if(clip.right < j) {
                        clip.right = j;
                    }
                    if(clip.top > i) {
                        clip.top = i;
                    }
                    if(clip.bottom < i) {
                        clip.bottom = i;
                    }
                }
            }
        }
        return result;
    }
}