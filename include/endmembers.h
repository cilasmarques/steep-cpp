#pragma once

#include "utils.h"
#include "constants.h"
#include "parameters.h"
#include "candidate.h"

/**
 * @brief  Compute the H0
 *
 * @param  net_radiation_line: Net radiation vector.
 * @param  soil_heat_flux: Soil heat flux vector.
 * @param  width_band: Band width.
 * @param  ho_line: H0 vector.
 *
 * @retval void
 */
void compute_H0(vector<double> net_radiation_line, vector<double> soil_heat_flux, int width_band, vector<double> &ho_line);

/**
 * @brief Calculates the four quartiles of a vector.
 *
 * @param target_vector: Vector to be calculated the quartiles.
 * @param v_quartile: Vector to store the quartiles.
 * @param height_band: Band height.
 * @param width_band: Band width.
 * @param first_interval: First interval.
 * @param middle_interval: Middle interval.
 * @param last_interval: Last interval.
 *
 * @retval void
 */
void get_quartiles(vector<vector<double>> target_vector, double *v_quartile, int height_band, int width_band, double first_interval, double middle_interval, double last_interval);

/**
 * @brief Check if the pixel is not a invalid value ( NaN or Inf )
 *
 * @param target_line: Line to be checked.
 * @param target_values: Vector to store the valid values.
 * @param width_band: Band width.
 * @param pos: Pointer to store the position of the valid values.
 */
void filter_valid_values(vector<double> target_line, double *target_values, int width_band, int *pos);

/**
 * @brief Get the hot and cold pixels based on the ASEBAL algorithm.
 *
 * @param ndvi_vector: NDVI vector.
 * @param surface_temperature_vector: Surface temperature vector.
 * @param albedo_vector: Albedo vector.
 * @param net_radiation_vector: Net radiation vector.
 * @param soil_heat_vector: Soil heat flux vector.
 * @param height_band: Band height.
 * @param width_band: Band width.
 *
 * @retval Candidate
 */
pair<Candidate, Candidate> getEndmembersASEBAL(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band);

/**
 * @brief Get the hot pixel based on the STEPP algorithm.
 *
 * @param ndvi_vector: NDVI vector.
 * @param surface_temperature_vector: Surface temperature vector.
 * @param albedo_vector: Albedo vector.
 * @param net_radiation_vector: Net radiation vector.
 * @param soil_heat_vector: Soil heat flux vector.
 * @param height_band: Band height.
 * @param width_band: Band width.
 *
 * @retval Candidate
 */
pair<Candidate, Candidate> getEndmembersSTEPP(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band);
