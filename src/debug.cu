#include "debug.h"

void printVector2x2(const vector<vector<double>> &matriz)
{
  for (const auto &linha : matriz)
  {
    for (const auto &elemento : linha)
    {
      std::cout << elemento << " ";
    }
    std::cout << std::endl;
  }
}

vector<vector<double>> tiffToVector(TIFF* tiffFile, int width, int height) {
    vector<vector<double>> data(height, vector<double>(width));

    for (int row = 0; row < height; ++row) {
        TIFFReadScanline(tiffFile, &data[row][0], row);
    }

    return data;
}
