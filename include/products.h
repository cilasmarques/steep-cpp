#pragma once

#include "utils.h"
#include "reader.h"
#include "candidate.h"
#include "constants.h"
#include "parameters.h"

/**
 * @brief  Struct to manage the products calculation.
 */
struct Products
{
  uint32_t width_band;
  uint32_t height_band;
  int nBytes_band;

  double H_pf_terra;
  double H_pq_terra;
  double rah_ini_pq_terra;
  double rah_ini_pf_terra;

  vector<vector<vector<double>>> radiance_vector;
  vector<vector<vector<double>>> reflectance_vector;
  vector<vector<double>> albedo_vector;
  vector<vector<double>> ndvi_vector;
  vector<vector<double>> soil_heat_vector;
  vector<vector<double>> surface_temperature_vector;
  vector<vector<double>> net_radiation_vector;
  vector<vector<double>> lai_vector;
  vector<vector<double>> evi_vector;
  vector<vector<double>> pai_vector;
  vector<vector<double>> enb_emissivity_vector;
  vector<vector<double>> eo_emissivity_vector;
  vector<vector<double>> ea_emissivity_vector;
  vector<vector<double>> short_wave_radiation_vector;
  vector<vector<double>> large_wave_radiation_surface_vector;
  vector<vector<double>> large_wave_radiation_atmosphere_vector;

  double *surface_temperature_pointer;
  double *d0_pointer;
  double *zom_pointer;
  double *ustar_pointer;
  double *kb1_pointer;
  double *aerodynamic_resistance_pointer;
  double *sensible_heat_flux_pointer;

  vector<vector<double>> sensible_heat_flux_vector;
  vector<vector<double>> latent_heat_flux_vector;
  vector<vector<double>> net_radiation_24h_vector;
  vector<vector<double>> evapotranspiration_fraction_vector;
  vector<vector<double>> sensible_heat_flux_24h_vector;
  vector<vector<double>> latent_heat_flux_24h_vector;
  vector<vector<double>> evapotranspiration_24h_vector;
  vector<vector<double>> evapotranspiration_vector;

  /**
   * @brief  Constructor.
   */
  Products();

  /**
   * @brief  Constructor.
   * @param  width_band: Band width.
   * @param  height_band: Band height.
   */
  Products(uint32_t width_band, uint32_t height_band);

  /**
   * @brief  The spectral radiance for each band is computed.
   * @param  bands_resampled: Resampled bands.
   * @param  width_band: Band width.
   * @param  sample_bands: Bands sample format.
   * @param  mtl: MTL struct.
   * @param  sensor: Sensor struct.
   * @param  line: Line to be computed.
   */
  void radiance_function(TIFF **bands_resampled, uint32_t width_band, uint16_t sample_bands, MTL mtl, Sensor sensor, int line);

  /**
   * @brief  The spectral reflectance for each band is computed.
   * @param  bands_resampled: Resampled bands.
   * @param  width_band: Band width.
   * @param  sample_bands: Bands sample format.
   * @param  mtl: MTL struct.
   * @param  sensor: Sensor struct.
   * @param  line: Line to be computed.
   */
  void reflectance_function(TIFF **bands_resampled, uint32_t width_band, uint16_t sample_bands, MTL mtl, Sensor sensor, int line);

  /**
   * @brief  The surface albedo is computed.
   * @param  tal_reader: TAL reader.
   * @param  sensor: Sensor struct.
   * @param  width_band: Band width.
   * @param  number_sensor: Number of sensor.
   * @param  line: Line to be computed.
   */
  void albedo_function(Reader tal_reader, Sensor sensor, uint32_t width_band, int number_sensor, int line);

  /**
   * @brief  The NDVI is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void ndvi_function(uint32_t width_band, int line);

  /**
   * @brief  The PAI is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void pai_function(uint32_t width_band, int line);

  /**
   * @brief  The LAI is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void lai_function(uint32_t width_band, int line);

  /**
   * @brief  The EVI is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void evi_function(uint32_t width_band, int line);

  /**
   * @brief  The emissivity is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void enb_emissivity_function(uint32_t width_band, int line);

  /**
   * @brief  The emissivity is computed.
   * @param  tal_reader: TAL reader.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void eo_emissivity_function(uint32_t width_band, int line);

  /**
   * @brief  The emissivity is computed.
   * @param  tal_reader: TAL reader.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void ea_emissivity_function(Reader tal_reader, uint32_t width_band, int line);

  /**
   * @brief  The surface temperature is computed.
   * @param  number_sensor: Number of sensor.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void surface_temperature_function(int number_sensor, uint32_t width_band, int line);

  /**
   * @brief  The short wave radiation is computed.
   * @param  tal_reader: TAL reader.
   * @param  mtl: MTL struct.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void short_wave_radiation_function(Reader tal_reader, MTL mtl, uint32_t width_band, int line);

  /**
   * @brief  The large wave radiation is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void large_wave_radiation_surface_function(uint32_t width_band, int line);

  /**
   * @brief  The large wave radiation is computed.
   * @param  width_band: Band width.
   * @param  temperature: Pixel's temperature.
   * @param  line: Line to be computed.
   */
  void large_wave_radiation_atmosphere_function(uint32_t width_band, double temperature, int line);

  /**
   * @brief  The net radiation is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void net_radiation_function(uint32_t width_band, int line);

  /**
   * @brief  The soil heat flux is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void soil_heat_flux_function(uint32_t width_band, int line);

  /**
   * @brief  The d0 is computed.
   * @param  line: Line to be computed.
   */
  void d0_fuction(int line);

  /**
   * @brief  The kb is computed.
   * @param  ndvi_max: Maximum NDVI.
   * @param  ndvi_min: Minimum NDVI.
   * @param  line: Line to be computed.
   */
  void kb_function(double ndvi_max, double ndvi_min, int line);

  /**
   * @brief  The zom is computed.
   * @param  A_ZOM: Coefficient A.
   * @param  B_ZOM: Coefficient B.
   * @param  line: Line to be computed.
   */
  void zom_fuction(double A_ZOM, double B_ZOM, int line);

  /**
   * @brief  The ustar is computed.
   * @param  u200: Wind speed at 200 m.
   * @param  zom_line: Zom line.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, int width_band, vector<double> &zom_line);

  /**
   * @brief  The ustar is computed.
   * @param  u10: Wind speed at 10 m.
   * @param  line: Line to be computed.
   */
  void ustar_fuction(double u10, int line);

  /**
   * @brief  The ustar is computed.
   * @param  u200: Wind speed at 200 m.
   * @param  zom_line: Zom line.
   * @param  width_band: Band width.
   * @param  ustar_line: Ustar line.
   */
  void ustar_fuction(double u200, vector<double> zom_line, int width_band, vector<double> &ustar_line);

  /**
   * @brief  The aerodynamic resistance is computed.
   * @param  line: Line to be computed.
   */
  void aerodynamic_resistance_fuction(int line);

  /**
   * @brief  The aerodynamic resistance is computed.
   * @param  ustar_line: Ustar line.
   * @param  width_band: Band width.
   * @param  aerodynamic_resistance_line: Aerodynamic resistance line.
   */
  void aerodynamic_resistance_fuction(vector<double> ustar_line, int width_band, vector<double> &aerodynamic_resistance_line);

  /**
   * @brief  The sensible heat flux is computed.
   * @param  a: Coefficient A.
   * @param  b: Coefficient B.
   * @param  line: Line to be computed.
   */
  void sensible_heat_flux_function(double a, double b, int line);

  /**
   * @brief  The latent heat flux is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void latent_heat_flux_function(int width_band, int line);

  /**
   * @brief  The net radiation is computed.
   * @param  Ra24h: Net radiation 24h.
   * @param  Rs24h: Solar radiation 24h.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void net_radiation_24h_function(double Ra24h, double Rs24h, int width_band, int line);

  /**
   * @brief  The evapotranspiration fraction is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void evapotranspiration_fraction_fuction(int width_band, int line);

  /**
   * @brief  The sensible heat flux is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void sensible_heat_flux_24h_fuction(int width_band, int line);

  /**
   * @brief  The latent heat flux is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void latent_heat_flux_24h_function(int width_band, int line);

  /**
   * @brief  The evapotranspiration is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void evapotranspiration_24h_function(Station station, int width_band, int line);

  /**
   * @brief  The evapotranspiration is computed.
   * @param  width_band: Band width.
   * @param  line: Line to be computed.
   */
  void evapotranspiration_function(int width_band, int line);

  /**
   * @brief  The  aerodynamic resistance convergence is computed.
   * @param  ndvi_min: Minimum NDVI.
   * @param  ndvi_max: Maximum NDVI.
   * @param  hot_pixel: Hot pixel.
   * @param  cold_pixel: Cold pixel.
   * @return  string: Time message.
   */
  string rah_correction_function_serial(double ndvi_min, double ndvi_max, Candidate hot_pixel, Candidate cold_pixel);

  /**
   * @brief  The  aerodynamic resistance convergence is computed.
   * @param  threads_num: Number of threads.
   * @param  ndvi_min: Minimum NDVI.
   * @param  ndvi_max: Maximum NDVI.
   * @param  hot_pixel: Hot pixel.
   * @param  cold_pixel: Cold pixel.
   * @return  string: Time message.
   */
  string rah_correction_function_threads(int threads_num, double ndvi_min, double ndvi_max, Candidate hot_pixel, Candidate cold_pixel);

  /**
   * @brief  The  aerodynamic resistance convergence is computed.
   * @param  ndvi_min: Minimum NDVI.
   * @param  ndvi_max: Maximum NDVI.
   * @param  hot_pixel: Hot pixel.
   * @param  cold_pixel: Cold pixel.
   * @return  string: Time message.
   */
  string rah_correction_function_blocks(double ndvi_min, double ndvi_max, Candidate hot_pixel, Candidate cold_pixel, int blocks_num);
};
