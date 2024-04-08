# pragma once

#include "parameters.h"
#include "candidate.h"

/**
 * @brief  Compute the initial value for the rah. (STEEP algorithm)
 * 
 * @param  station: Station struct.
 * @param  start_line: Start line.
 * @param  end_line: End line.
 * @param  ndvi_min: Minimum NDVI value.
 * @param  ndvi_max: Maximum NDVI value.
 * @param  u10: Wind speed at 10 m.
*/
void rah_initial_value_STEEP(Station station, int start_line, int end_line, double ndvi_min, double ndvi_max, double u10);

/**
 * @brief  Compute the rah correction cycle. (STEEP algorithm)
 * 
 * @param  start_line: Start line.
 * @param  end_line: End line.
 * @param  hot_pixel: Hot pixel.
 * @param  cold_pixel: Cold pixel.
 * @param  a: Coefficient A.
 * @param  b: Coefficient B.
*/
void rah_correction_cycle_STEEP(int start_line, int end_line, Candidate hot_pixel, Candidate cold_pixel, double a, double b);

/**
 * @brief  Compute the final value for the rah. 
 * 
 * @param  start_line: Start line.
 * @param  end_line: End line.
 * @param  a: Coefficient A.
 * @param  b: Coefficient B.
*/
void sensible_heat_flux_final(int start_line, int end_line, double a, double b);
