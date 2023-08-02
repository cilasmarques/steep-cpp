#include "utils.h"

void printVector2x2(const std::vector<std::vector<double>> &matriz);
void printTIF2x2(TIFF *tiff, uint32 height_band, uint32 width_band);
void vectorToArray(vector<vector<double>> vec, double matriz[][35]);
void tiffToArray(const char* filename, double** matriz, uint32& imageWidth, uint32& imageLength);
std::vector<std::vector<double>> tiffToVector(TIFF* tiffFile, int width, int height);
std::vector<std::vector<double>> arrayToVector(double matriz[][35], int height_band, int width_band);
