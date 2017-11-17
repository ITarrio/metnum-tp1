#include <iostream>
#include <fstream>
#include "utils.h"
#include "reconstruct3d.h"

#define SAVE_NORMAL true

#define DEFAULT_OPTS options(false)

std::vector<matrix<double>> loadGrayImages(const std::vector<std::string>& files) {
    std::vector<matrix<double>> res;
    for (auto it = files.begin(); it != files.end(); ++it) {
        res.push_back(Utils::loadGrayImage(*it));
    }
    return res;
}

void printUsage() {
    std::cerr << "Uso: ./procesar <lista_referencia.txt> <lista_figura.txt>\n"
              << "\n"
              << "Ejemplo: ./procesar mate.txt caballo.txt\n"
              << "\n"
              << "(el archivo lista_referencia.txt debe listar\n"
              << "las imágenes de la esfera para usar al calibrar)\n"
              << std::endl;
}

std::string getSubfolder(const std::string& path) {
    size_t pos = 0;
    while(path.find('/', pos) != std::string::npos) {
        pos = path.find('/', pos) + 1;
    }
    return path.substr(0, pos);
}

std::vector<std::string> getFileNames(const std::string &file, std::string &mask) {
    std::vector<std::string> names;
    std::ifstream list(file);

    std::string folder = getSubfolder(file);
    uint n;
    list >> n;
    for (uint i = 0; i < n; ++i) {
        std::string name;
        list >> name;
        names.push_back(folder + name);
    }
    list >> mask;
    mask = folder + mask;
    return names;
}

std::string getOutputName(const std::string &file) {
    size_t pos = 0;
    while(file.find('/', pos) != std::string::npos) {
        pos = file.find('/', pos) + 1;
    }
    size_t ext = pos;
    while(file.find('.', ext + 1) != std::string::npos) {
        ext = file.find('.', ext + 1);
    }
    return file.substr(pos, ext - pos);
}

int main(int argc, char** argv) {
    if(argc != 3) {
        printUsage();
        return 1;
    }
    std::string referenceName(argv[1]);
    std::string figureName(argv[2]);
    //----CALIBRATION
    std::string referenceMaskName;
    std::vector<std::string> referenceFiles = getFileNames(referenceName, referenceMaskName);
    matrix<double> referenceMask = Utils::loadGrayImage(referenceMaskName);
    vector<matrix<double>> references = loadGrayImages(referenceFiles);
    for (matrix<double> ref : references) {
        ref = Mask::apply_mask(ref, referenceMask);
    }
    vector<direction> calibrations = Calibration::calibrate(references, referenceMask);
    //----LOADING
    std::string figureMaskName;
    std::vector<std::string> figureFiles = getFileNames(figureName, figureMaskName);
    matrix<double> mask = Utils::loadGrayImage(figureMaskName);
    vector<matrix<double>> imgs = loadGrayImages(figureFiles);
    for (matrix<double> img : references) {
        img = Mask::apply_mask(img, mask);
    }
    //----NORMAL-FIELD
    // TODO: CHOOSE THE BEST 3
    std::string outputName = getOutputName(figureName);
    matrix<row<double>> normal = normalField(
            imgs[0], imgs[1], imgs[2],
            calibrations[0], calibrations[1], calibrations[2], DEFAULT_OPTS);
#if SAVE_NORMAL
        Utils::saveMatrix3dFiles(normal, outputName);
#endif
    //----DEPTH
    matrix<double> depth = findDepth(normal);
    //Utils::saveMatrixFile(depth, outputName + ".depth.csv");
    //----END
    return 0;
}
