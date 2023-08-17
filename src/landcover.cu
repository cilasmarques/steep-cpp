#include "landcover.h"

bool checkLandCode(int value)
{
  return (value == AGP) || (value == PAS) || (value == AGR) || (value == CAP) || (value == CSP) || (value == MAP);
}

void testLandCoverHomogeneity(TIFF *landCover, vector<vector<double>> &mask_lc_vector)
{
  uint16 sample_bands;
  uint32 height_band, width_band;
  TIFFGetField(landCover, TIFFTAG_IMAGELENGTH, &height_band);
  TIFFGetField(landCover, TIFFTAG_IMAGEWIDTH, &width_band);
  TIFFGetField(landCover, TIFFTAG_SAMPLEFORMAT, &sample_bands);

  double **buffer = (double **)malloc(7 * sizeof(double *));

  for (int i = 0; i < 7; i++)
  {
    buffer[i] = (double *)malloc(width_band * sizeof(double));
  }

  int relation[7] = {-1, -1, -1, -1, -1, -1, -1}, aux;

  for (int line = 0; line < height_band; line++)
  {
    for (int column = 0; column < width_band; column++)
    {
      int pixel_value;
      aux = line % 7;

      if (relation[aux] != line)
      {
        TIFFReadScanline(landCover, buffer[aux], line);
        relation[aux] = line;
      }

      mask_lc_vector[line][column] = false;
      pixel_value = buffer[aux][column];
      if (checkLandCode(pixel_value))
      { // Verify if the pixel is an AGR pixel
        mask_lc_vector[line][column] = true;

        for (int i = -3; i <= 3 && mask_lc_vector[line][column]; i++)
        {
          for (int j = -3; j <= 3 && mask_lc_vector[line][column]; j++)
          {
            // Check if the neighbor is AGR too
            if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band)
            {
              aux = (line + j) % 7;
              if (relation[aux] != (line + j))
              {
                TIFFReadScanline(landCover, buffer[aux], line + j);
                relation[aux] = (line + j);
              }

              pixel_value = buffer[aux][column + i];
              if (!isnan(pixel_value))
                if (!checkLandCode(pixel_value))
                  mask_lc_vector[line][column] = false;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < 7; i++)
  {
    free(buffer[i]);
  }
  free(buffer);
}

void testHomogeneity(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> mask_lc_vector, int height_band, int width_band, vector<vector<double>> &homo_vector)
{
  int relation[7] = {-1, -1, -1, -1, -1, -1, -1}, aux;

  for (int line = 0; line < height_band; line++)
  {
    for (int column = 0; column < width_band; column++)
    {
      if (mask_lc_vector[line][column] == true)
      { // Verify if the pixel passed the land cover test

        aux = line % 7;
        double pixel_value;

        vector<double> ndvi_neighbors;
        vector<double> ts_neighbors;
        vector<double> albedo_neighbors;

        if (relation[aux] != line)
          relation[aux] = line;

        if (!isnan(ndvi_vector[aux][column]))
        {
          for (int i = -3; i <= 3; i++)
          {
            for (int j = -3; j <= 3; j++)
            {
              // Add for the NDVI, TS and Albedo the value of neighbors pixels into the respective vector
              if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band)
              {

                aux = (line + j) % 7;

                if (relation[aux] != (line + j))
                  relation[aux] = line + j;

                pixel_value = ndvi_vector[aux][column + i];
                if (!isnan(pixel_value))
                  ndvi_neighbors.push_back(pixel_value);

                pixel_value = surface_temperature_vector[aux][column + i];
                if (!isnan(pixel_value))
                  ts_neighbors.push_back(pixel_value);

                pixel_value = albedo_vector[aux][column + i];
                if (!isnan(pixel_value))
                  albedo_neighbors.push_back(pixel_value);
              }
            }
          }

          // Do the calculation of the dispersion measures from the NDVI, TS and Albedo
          double meanNDVI, meanTS, meanAlb;
          double sdNDVI, sdTS, sdAlb;
          double cvNDVI, cvAlb;
          double sumNDVI = 0, sumTS = 0, sumAlb = 0;

          for (int i = 0; i < ndvi_neighbors.size(); i++)
          {
            sumNDVI += ndvi_neighbors[i];
            sumTS += ts_neighbors[i];
            sumAlb += albedo_neighbors[i];
          }

          meanNDVI = sumNDVI / ndvi_neighbors.size();
          meanTS = sumTS / ts_neighbors.size();
          meanAlb = sumAlb / albedo_neighbors.size();

          sumNDVI = 0, sumTS = 0, sumAlb = 0;

          for (int i = 0; i < ndvi_neighbors.size(); i++)
          {
            sumNDVI += (ndvi_neighbors[i] - meanNDVI) * (ndvi_neighbors[i] - meanNDVI);
            sumTS += (ts_neighbors[i] - meanTS) * (ts_neighbors[i] - meanTS);
            sumAlb += (albedo_neighbors[i] - meanAlb) * (albedo_neighbors[i] - meanAlb);
          }

          sdNDVI = sqrt(sumNDVI / ndvi_neighbors.size());
          sdTS = sqrt(sumTS / ts_neighbors.size());
          sdAlb = sqrt(sumAlb / albedo_neighbors.size());

          cvNDVI = sdNDVI / meanNDVI;
          cvAlb = sdAlb / meanAlb;

          // Check if the pixel is eligible
          homo_vector[line][column] = (cvNDVI < 0.25) && (cvAlb < 0.25) && (sdTS < 1.5);
        }
        else
        {
          homo_vector[line][column] = false;
        }
      }
    }
  }
}

void testMorphological(vector<vector<double>> homo_vector, int groupSize, int height_band, int width_band, vector<vector<double>> &morph_vector)
{
  // Apply the routine
  queue<pair<int, int>> fila;
  set<pair<int, int>> cont;

  for (int line = 0; line < height_band; line++)
  {
    for (int col = 0; col < width_band; col++)
    {
      if (homo_vector[line][col] == 1)
      {
        fila.push({line, col});
        cont.insert({line, col});
        homo_vector[line][col] = -1;

        while (!fila.empty())
        {
          int i = fila.front().first;
          int j = fila.front().second;
          fila.pop();

          if (j + 1 < width_band)
          {
            if (homo_vector[i][j + 1] == 1)
            {
              fila.push({i, j + 1});
              cont.insert({i, j + 1});
              homo_vector[i][j + 1] = -1;
            }

            if (i + 1 < height_band && homo_vector[i + 1][j + 1] == 1)
            {
              fila.push({i + 1, j + 1});
              cont.insert({i + 1, j + 1});
              homo_vector[i + 1][j + 1] = -1;
            }

            if (i > 0 && homo_vector[i - 1][j + 1] == 1)
            {
              fila.push({i - 1, j + 1});
              cont.insert({i - 1, j + 1});
              homo_vector[i - 1][j + 1] = -1;
            }
          }

          if (j > 0)
          {
            if (homo_vector[i][j - 1] == 1)
            {
              fila.push({i, j - 1});
              cont.insert({i, j - 1});
              homo_vector[i][j - 1] = -1;
            }

            if (i + 1 < height_band && homo_vector[i + 1][j - 1] == 1)
            {
              fila.push({i + 1, j - 1});
              cont.insert({i + 1, j - 1});
              homo_vector[i + 1][j - 1] = -1;
            }

            if (i > 0 && homo_vector[i - 1][j - 1] == 1)
            {
              fila.push({i - 1, j - 1});
              cont.insert({i - 1, j - 1});
              homo_vector[i - 1][j - 1] = -1;
            }
          }

          if (i + 1 < height_band && homo_vector[i + 1][j] == 1)
          {
            fila.push({i + 1, j});
            cont.insert({i + 1, j});
            homo_vector[i + 1][j] = -1;
          }

          if (i > 0 && homo_vector[i - 1][j] == 1)
          {
            fila.push({i - 1, j});
            cont.insert({i - 1, j});
            homo_vector[i - 1][j] = -1;
          }
        }

        int group = cont.size();

        for (auto elem : cont)
        {
          morph_vector[elem.first][elem.second] = (group >= groupSize);
        }

        cont.clear();
      }
      else if (homo_vector[line][col] == 0)
      {
        morph_vector[line][col] = 0;
      }
    }
  }
}
