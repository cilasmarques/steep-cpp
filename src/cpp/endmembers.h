#pragma once

#include "utils.h"
#include "reader.h"
#include "landcover.h"
#include "constants.h"
#include "parameters.h"
#include "candidate.h"

void compute_H0(vector<double> net_radiation_line, vector<double> soil_heat_flux, int width_band, vector<double> &ho_line);
void get_quartiles(vector<vector<double>> target_vector, double *v_quartile, int height_band, int width_band, double first_interval, double mid_interval, double last_interval);
void filter_valid_values(vector<double> target_line, double *target_values, int width_band, int *pos);

pair<Candidate, Candidate> getColdHotPixelsSTEPP(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band);
