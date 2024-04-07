#pragma once

#include "constants.h"
#include "utils.h"
#include "reader.h"
#include "parameters.h"
#include "candidate.h"

bool checkLandCode(int value);
void testLandCoverHomogeneity(TIFF *landCover, vector<vector<double>> &mask_vector);
void testHomogeneity(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> maskLC, int height_band, int width_band, vector<vector<double>> &outputH);
void testMorphological(vector<vector<double>> homo_vector, int groupSize, int height_band, int width_band, vector<vector<double>> &morph_vector);
