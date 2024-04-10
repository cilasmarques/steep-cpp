#pragma once

#include "constants.h"

/**
 * @brief  Compute the rah correction cycle. (STEEP algorithm)
 * 
 * @param surface_temperature_pointer  Surface temperature
 * @param d0_pointer  Zero plane displacement height
 * @param kb1_pointer  KB-1 stability parameter
 * @param zom_pointer  Roughness length for momentum
 * @param ustarR_pointer  Ustar pointer for reading
 * @param ustarW_pointer  Ustar pointer for writing
 * @param rahR_pointer  Rah pointer for reading
 * @param rahWL_pointer  Rah pointer for writing
 * @param H_pointer  Sensible heat flux
 * @param a  Coefficient a
 * @param b  Coefficient b
 * @param height  Height of the input data
 * @param width  Width of the input data
*/
__global__ void rah_correction_cycle_STEEP(double *surface_temperature_pointer, double *d0_pointer, double *kb1_pointer, double *zom_pointer, double *ustarR_pointer,
                                           double *ustarW_pointer, double *rahR_pointer, double *rahWL_pointer, double *H_pointer, double a, double b, int height,
                                           int width);
