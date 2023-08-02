#include "debug.h"

void printVector2x2(const std::vector<std::vector<double>> &matriz)
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

void printTIF2x2(TIFF *tiff, uint32 height_band, uint32 width_band)
{
  double tiff_line[width_band];
  for (int line = 0; line < height_band; ++line)
  {
    read_line_tiff(tiff, tiff_line, line); 
    for (int col = 0; col < width_band; col++) {
      cout << tiff_line[col] << " ";
    }
    cout << endl;
  }
}

void vectorToArray(vector<vector<double>> vec, double matriz[][35]) {
    for (int i = 0; i < 35; i++) {
        for (int j = 0; j < 35; j++) {
            matriz[i][j] = vec[i][j];
        }
    }
}

std::vector<std::vector<double>> arrayToVector(double matriz[][35], int height_band, int width_band) {
    std::vector<std::vector<double>> vetorDeVetores;

    for (int i = 0; i < height_band; i++) {
        std::vector<double> vetorInterno;
        for (int j = 0; j < width_band; j++) {
            vetorInterno.push_back(matriz[i][j]);
        }
        vetorDeVetores.push_back(vetorInterno);
    }

    return vetorDeVetores;
}

std::vector<std::vector<double>> tiffToVector(TIFF* tiffFile, int width, int height) {
    std::vector<std::vector<double>> data(height, std::vector<double>(width));

    for (int row = 0; row < height; ++row) {
        TIFFReadScanline(tiffFile, &data[row][0], row);
    }

    return data;
}

void tiffToArray(const char* filename, double** matriz, uint32& imageWidth, uint32& imageLength) {
    TIFF* tif = TIFFOpen(filename, "r");
    uint32* raster;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
    size_t npixels = imageWidth * imageLength;
    raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));
    TIFFReadRGBAImage(tif, imageWidth, imageLength, raster, 0);

    for (uint32 row = 0; row < imageLength; row++) {
        for (uint32 col = 0; col < imageWidth; col++) {
            matriz[row][col] = (double)raster[(imageLength - row - 1) * imageWidth + col];
        }
    }

    _TIFFfree(raster);
    TIFFClose(tif);
}