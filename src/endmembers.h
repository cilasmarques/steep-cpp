#pragma once

#include "utils.h"
#include "reader.h"
#include "landcover.h"
#include "constants.h"
#include "parameters.h"
#include "candidate.h"

void compute_H0(vector<double> net_radiation_line, vector<double> soil_heat_flux, int width_band, vector<double> &ho_line);
void get_quartiles(vector<vector<double>> target_matriz, double *v_quartile, int height_band, int width_band, double first_interval, double last_interval);
void filter_valid_values(vector<double> target_line, double *target_values, int width_band, int *pos);

Candidate getHotPixelASEBAL(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band);
Candidate getColdPixelASEBAL(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band);

Candidate getHotPixelSTEPP(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band);
Candidate getColdPixelSTEPP(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band);

pair<Candidate, Candidate> getColdHotPixelsESA(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, TIFF** landCover, int height_band, int width_band, string output_path);