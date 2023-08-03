#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

void compute_H0(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);

void filter_valid_values(double *target_line, double *target_values, int width_band, int *pos);

void get_quartiles(std::vector<std::vector<double>> target_matriz, double *v_quartile, int height_band, int width_band, double first_interval, double last_interval);

Candidate getHotPixelSTEPP(std::vector<std::vector<double>> ndvi_matriz, std::vector<std::vector<double>> surface_temperature_matriz, std::vector<std::vector<double>> albedo_matriz, std::vector<std::vector<double>> net_radiation_matriz, std::vector<std::vector<double>> soil_heat_matriz, int height_band, int width_band);

Candidate getColdPixelSTEPP(std::vector<std::vector<double>> ndvi_matriz, std::vector<std::vector<double>> surface_temperature_matriz, std::vector<std::vector<double>> albedo_matriz, std::vector<std::vector<double>> net_radiation_matriz, std::vector<std::vector<double>> soil_heat_matriz, int height_band, int width_band);
