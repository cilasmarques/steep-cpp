#pragma once

#include "constants.h"

/**
 * @brief  Compute the rah correction cycle. (STEEP algorithm)
 *
 * @param  start_line: The start line of the band.
 * @param  end_line: The end line of the band.
 * @param  width_band: The width of the band.
 * @param  a: The a parameter.
 * @param  b: The b parameter.
 * @param  surface_temperature_pointer: The surface temperature vector.
 * @param  d0_pointer: The d0 vector.
 * @param  aerodynamic_resistance_previous: The aerodynamic resistance vector of the previous iteration.
 * @param  ustar_previous: The ustar vector of the previous iteration.
 * @param  zom_pointer: The zom vector.
 * @param  kb1_pointer: The kb1 vector.
 * @param  sensible_heat_flux_pointer: The sensible heat flux vector.
 * @param  ustar_pointer: The ustar vector.
 * @param  aerodynamic_resistance_pointer: The aerodynamic resistance vector.
 */
void rah_correction_cycle_STEEP(int start_line, int end_line, int width_band, double a, double b, double *surface_temperature_pointer,
                                double *d0_pointer, double *zom_pointer, double *kb1_pointer, double *sensible_heat_flux_pointer, 
                                double *ustar_pointer, double *aerodynamic_resistance_pointer);
