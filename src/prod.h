#pragma once

#include "types.h"
#include "utils.h"
#include "tiffs.h"
#include "parameters.h"
#include "pixel_reader.h"
#include <fstream>
#include <iostream>

/**
 * @brief  The spectral radiance for each band is computed.
 * @param  bands_resampled[]: Satellite bands.
 * @param  mtl: MTL struct.
 * @param  sensor: Sensor struct.
 * @param  width_band: Band width.
 * @param  line: Line to be calculated.
 * @param  radiance_line[][8]: Auxiliary array for save the calculated value of radiance for each band.
 */
void radiance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, double radiance_line[][8]);
void reflectance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, double reflectance_line[][8]);
void albedo_function(PixelReader tal_reader, Sensor sensor, uint32 width_band, int number_sensor, double reflectance_line[][8], double albedo_line[]);
void ndvi_function(double reflectance_line[][8], uint32 width_band, double ndvi_line[]);

/**
 * @brief  Plant Area Index (PAI) is computed.
 * @note   Values for PAI range between -1 and 1.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  pai_line[]: Auxiliary array for save the calculated value of NDVI for the line.
 */
void pai_function(double reflectance_line[][8], uint32 width_band, double pai_line[]);

/**
 * @brief  Leaf Area Index (LAI) is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  lai_line[]: Auxiliary array for save the calculated value of LAI for the line.
 */
void lai_function(double reflectance_line[][8], uint32 width_band, double lai_line[]);

/**
 * @brief  Enhanced Vegetation Index (EVI) is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  evi_line[]: Auxiliary array for save the calculated value of EVI for the line.
 */
void evi_function(double reflectance_line[][8], uint32 width_band, double evi_line[]);

/**
 * @brief  Calculates emissivity representing surface behavior for thermal emission in the relatively narrow band 6 of Landsat (10.4 to 12.5 µm),
           expressed as enb.
 * @param  lai_line[]: Array containing the specified line from the LAI computation.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  enb_emissivity_line[]: Auxiliary array for save the calculated value of Enb for the line.
 */
void enb_emissivity_function(double lai_line[], double ndvi_line[], uint32 width_band, double enb_emissivity_line[]);

/**
 * @brief  Calculates emissivity representing surface behavior for thermal emission in the broad thermal spectrum (6 to 14 µm), expressed as eο.  
 * @param  lai_line[]: Array containing the specified line from the LAI computation.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  eo_emissivity_line[]: Auxiliary array for save the calculated value of Eo for the line.
 */
void eo_emissivity_function(double lai_line[], double ndvi_line[], uint32 width_band, double eo_emissivity_line[]);

/**
 * @brief  Calculates the atmospheric emissivity (ea).
 * @param  tal_reader: Object containing the specified line from the tal computation.
 * @param  width_band: Band width.
 * @param  ea_emissivity_line[]: Auxiliary array for save the calculated value of Ea for the line.
 */
void ea_emissivity_function(PixelReader tal_reader, uint32 width_band, double ea_emissivity_line[]);

/**
 * @brief  The surface temperature (TS) is computed.
 * @param  radiance_line[][8]: Radiance for the specific line for each band.
 * @param  enb_emissivity_line[]: Array containing the specified line from the Enb computation.
 * @param  number_sensor: Number of the satellite sensor.
 * @param  width_band: Band width.
 * @param  surface_temperature_line[]: Auxiliary array for save the calculated value of TS for the line.
 */
void surface_temperature_function(double radiance_line[][8], double enb_emissivity_line[], int number_sensor, uint32 width_band, double surface_temperature_line[]);

/**
 * @brief  Computes Short Wave Radiation (Rs).
 * @param  tal_reader: Object containing the specified line from the tal computation.
 * @param  mtl: MTL Struct.
 * @param  width_band: Band width.
 * @param  short_wave_radiation_line[]: Auxiliary array for save the calculated value of Rs for the line.
 */
void short_wave_radiation_function(PixelReader tal_reader, MTL mtl, uint32 width_band, double short_wave_radiation_line[]);

/**
 * @brief  Computes Large Wave Radiation from Surface (RLSup)
 * @param  eo_emissivity_line[]: Array containing the specified line from the Eo computation.
 * @param  surface_temperature_line[]: Array containing the specified line from the TS computation.
 * @param  width_band: Band width.
 * @param  large_wave_radiation_surface_line[]: Auxiliary array for save the calculated value of RLSup for the line.
 */
void large_wave_radiation_surface_function(double eo_emissivity_line[], double surface_temperature_line[], uint32 width_band, double large_wave_radiation_surface_line[]);

/**
 * @brief  Computes Large Wave Radiation from Atmosphere (RLatm)
 * @param  ea_emissivity_line[]: Array containing the specified line from the Ea computation.
 * @param  width_band: Band width.
 * @param  temperature: Near surface air temperature in Kelvin.
 * @param  large_wave_radiation_atmosphere_line[]: Auxiliary array for save the calculated value of RLatm for the line.
 */
void large_wave_radiation_atmosphere_function(double ea_emissivity_line[], uint32 width_band, double temperature, double large_wave_radiation_atmosphere_line[]);

/**
 * @brief  The net surface radiation flux (Rn) is computed.
 * @param  short_wave_radiation_line[]: Array containing the specified line from the Rs computation.
 * @param  large_wave_radiation_surface_line[]: Array containing the specified line from the RLSup computation.
 * @param  large_wave_radiation_atmosphere_line[]: Array containing the specified line from the RLatm computation.
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  eo_emissivity_line[]: Array containing the specified line from the Eo computation.
 * @param  width_band: Band width.
 * @param  net_radiation_line[]: Auxiliary array for save the calculated value of Rn for the line.
 */
void net_radiation_function(double short_wave_radiation_line[], double large_wave_radiation_surface_line[],
                            double large_wave_radiation_atmosphere_line[], double albedo_line[],
                            double eo_emissivity_line[], uint32 width_band, double net_radiation_line[]);

/**
 * @brief  Computes the Soil heat flux (G).    
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  surface_temperature_line[]: Array containing the specified line from the TS computation.
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  width_band: Band width.
 * @param  soil_heat_flux[]: Auxiliary array for save the calculated value of G for the line.
 */
void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], uint32 width_band, double soil_heat_flux[]);
