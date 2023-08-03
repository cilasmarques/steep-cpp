#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

void d0_fuction(std::vector<double> pai_line, int width_band, double d0_line[]);
void kb_function(double ustar_line[], double zom_line[], std::vector<double> pai_line, double ndvi_line[], double ndvi_max, double ndvi_min, int width_band, double kb1_line[]);

/**
 * @brief  Computes the momentum roughness length (zom).
 * @param  A_ZOM: Correlation constant a.
 * @param  B_ZOM: Correlation constant b.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  zom_line[]: Auxiliary array for save the calculated value of zom for the line.
 */
void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[], double d0_line[], std::vector<double> pai_line);

/**
 * @brief  The friction velocity (u*) is computed.
 * @param  u10: Wind speed at 200 m.
 * @param  zom_line[]: Array containing the specified line from the zom computation.
 * @param  width_band: Band width.
 * @param  ustar_line[]: Auxiliary array for save the calculated value of ustar for the line.
 */
void ustar_fuction(double u10, double zom_line[], int width_band, double ustar_line[], double d0_line[]);

/**
 * @brief  Computes the aerodynamic resistance (Rah).   
 * @param  ustar_line[]: Array containing the specified line from the ustar computation.
 * @param  width_band: Band width.
 * @param  aerodynamic_resistance_line[]: Auxiliary array for save the calculated value of Rah for the line.
 */
void aerodynamic_resistance_fuction(double ustar_line[], int width_band, double aerodynamic_resistance_line[], double zom_line[], double d0_line[], double kb1_line[]);

/**
 * @brief  Computes Latent Heat Flux (LE).  
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_flux_line[]: Array containing the specified line from the G computation.
 * @param  sensible_heat_flux_line[]: Array containing the specified line from the H computation.
 * @param  width_band: Band width.
 * @param  latent_heat_flux[]: Auxiliary array for save the calculated value of LE for the line.
 */
void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[]);

/**
 * @brief  Calculates the Net Radiation for 24 hours (Rn24h).
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  Ra24h: Ra24h: Extraterrestrial Radiation defined as solar short wave radiation in the absence of an atmosphere (Ra24h).
 * @param  Rs24h: Short wave radiation incident in 24 hours (Rs24h).
 * @param  width_band: Band width.
 * @param  net_radiation_24h_line[]: Auxiliary array for save the calculated value of Rn24h for the line.
 */
void net_radiation_24h_function(double albedo_line[], double Ra24h, double Rs24h, int width_band, double net_radiation_24h_line[]);

/**
 * @brief  The Reference ET Fraction (EF) is computed.
 * @param  latent_heat_flux_line[]: Array containing the specified line from the LE computation.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_line[]: Array containing the specified line from the G computation.
 * @param  width_band: Band width.
 * @param  evapotranspiration_fraction_line[]: Auxiliary array for save the calculated value of EF for the line.
 */
void evapotranspiration_fraction_fuction(double latent_heat_flux_line[], double net_radiation_line[], double soil_heat_line[], int width_band, double evapotranspiration_fraction_line[]);

/**
 * @brief  Computes Sensible Heat Flux for 24 hours (H24h).
 * @param  evapotranspiration_fraction_line[]: Array containing the specified line from the EF computation.
 * @param  net_radiation_24h_line[]: Array containing the specified line from the Rn24h computation.
 * @param  width_band: Band width.
 * @param  sensible_heat_flux_24h_line[]: Auxiliary array for save the calculated value of H24h for the line.
 */
void sensible_heat_flux_24h_fuction(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double sensible_heat_flux_24h_line[]);

/**
 * @brief  Calculates Latente Heat Flux for 24 hours (LE24h).
 * @param  evapotranspiration_fraction_line[]: Array containing the specified line from the EF computation.
 * @param  net_radiation_24h_line[]: Array containing the specified line from the Rn24h computation.
 * @param  width_band: Band width.
 * @param  latent_heat_flux_24h_line[]: Auxiliary array for save the calculated value of LE24h for the line.
 */
void latent_heat_flux_24h_function(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double latent_heat_flux_24h_line[]);

/**
 * @brief  Computes the Evapotranspiration for 24 hours (ET24h)
 * @param  latent_heat_flux_24h_line[]: Array containing the specified line from the LE24h computation.
 * @param  station: Station struct.
 * @param  width_band: Band width.
 * @param  evapotranspiration_24h_line[]: Auxiliary array for save the calculated value of ET24h for the line.
 */
void evapotranspiration_24h_function(double latent_heat_flux_24h_line[], Station station, int width_band, double evapotranspiration_24h_line[]);
