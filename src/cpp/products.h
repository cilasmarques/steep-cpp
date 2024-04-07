#pragma once

#include "debug.h"
#include "constants.h"
#include "utils.h"
#include "reader.h"
#include "candidate.h"
#include "parameters.h"

/**
 * @brief  Struct to manage the products calculation.
 */
struct Products
{
  uint32_t width_band;
  uint32_t height_band;

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
  vector<vector<double>> d0_vector;
  vector<vector<double>> zom_vector;
  vector<vector<double>> ustar_vector;
  vector<vector<double>> kb1_vector;
  vector<vector<double>> aerodynamic_resistance_vector;
  vector<vector<double>> sensible_heat_flux_vector;
  vector<vector<double>> ustar_previous;
  vector<vector<double>> aerodynamic_resistance_previous;
  vector<vector<double>> latent_heat_flux_vector;
  vector<vector<double>> net_radiation_24h_vector;
  vector<vector<double>> evapotranspiration_fraction_vector;
  vector<vector<double>> sensible_heat_flux_24h_vector;
  vector<vector<double>> latent_heat_flux_24h_vector;
  vector<vector<double>> evapotranspiration_24h_vector;
  vector<vector<double>> evapotranspiration_vector;

  Products(uint32_t width_band, uint32_t height_band);

  void radiance_function(TIFF **bands_resampled, uint32_t width_band, uint16_t sample_bands, MTL mtl, Sensor sensor, int line);
  void reflectance_function(TIFF **bands_resampled, uint32_t width_band, uint16_t sample_bands, MTL mtl, Sensor sensor, int line);
  void albedo_function(Reader tal_reader, Sensor sensor, uint32_t width_band, int number_sensor, int line);
  void ndvi_function(uint32_t width_band, int line);
  void pai_function(uint32_t width_band, int line);
  void lai_function(uint32_t width_band, int line);
  void evi_function(uint32_t width_band, int line);
  void enb_emissivity_function(uint32_t width_band, int line);
  void eo_emissivity_function(uint32_t width_band, int line);
  void ea_emissivity_function(Reader tal_reader, uint32_t width_band, int line);
  void surface_temperature_function(int number_sensor, uint32_t width_band, int line);
  void short_wave_radiation_function(Reader tal_reader, MTL mtl, uint32_t width_band, int line);
  void large_wave_radiation_surface_function(uint32_t width_band, int line);
  void large_wave_radiation_atmosphere_function(uint32_t width_band, double temperature, int line);
  void net_radiation_function(uint32_t width_band, int line);
  void soil_heat_flux_function(uint32_t width_band, int line);

  void d0_fuction(int line);
  void kb_function(double ndvi_max, double ndvi_min, int line);
  void zom_fuction(double A_ZOM, double B_ZOM, int line);
  void zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, int width_band, vector<double> &zom_line);
  void ustar_fuction(double u10, int line);
  void ustar_fuction(double u200, vector<double> zom_line, int width_band, vector<double> &ustar_line);
  void aerodynamic_resistance_fuction(vector<double> ustar_line, int width_band, vector<double> &aerodynamic_resistance_line);
  void aerodynamic_resistance_fuction(int line);

  void sensible_heat_function_STEEP(Candidate hot_pixel, Candidate cold_pixel, Station station, uint32_t height_band, uint32_t width_band, int threads_num);
  void sensible_heat_function_default(Candidate hot_pixel, Candidate cold_pixel, Station station, uint32_t height_band, uint32_t width_band, int threads_num);

  void latent_heat_flux_function(int width_band, int line);
  void net_radiation_24h_function(double Ra24h, double Rs24h, int width_band, int line);
  void evapotranspiration_fraction_fuction(int width_band, int line);
  void sensible_heat_flux_24h_fuction(int width_band, int line);
  void latent_heat_flux_24h_function(int width_band, int line);
  void evapotranspiration_24h_function(Station station, int width_band, int line);
  void evapotranspiration_function(int width_band, int line);

  void rah_correction_cycle_STEEP(int start_line, int end_line, Candidate hot_pixel, Candidate cold_pixel, double a, double b);
  void rah_initial_value_STEEP(Station station, int start_line, int end_line, double ndvi_min, double ndvi_max, double u10);
  void sensible_heat_flux_final(int start_line, int end_line, double a, double b);
};
