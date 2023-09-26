#pragma once

#include "debug.h"
#include "constants.h"
#include "utils.h"
#include "reader.h"
#include "candidate.h"
#include "parameters.h"

void radiance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, vector<vector<double>> &radiance_line);
void reflectance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, vector<vector<double>> &reflectance_line);
void albedo_function(Reader tal_reader, Sensor sensor, uint32 width_band, int number_sensor, vector<vector<double>> reflectance_line, vector<double> &albedo_line);
void ndvi_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &ndvi_line);
void pai_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &pai_line);
void lai_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &lai_line);
void evi_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &evi_line);
void enb_emissivity_function(vector<double> lai_line, vector<double> ndvi_line, uint32 width_band, vector<double> &enb_emissivity_line);
void eo_emissivity_function(vector<double> lai_line, vector<double> ndvi_line, uint32 width_band, vector<double> &eo_emissivity_line);
void ea_emissivity_function(Reader tal_reader, uint32 width_band, vector<double> &ea_emissivity_line);
void surface_temperature_function(vector<vector<double>> radiance_line, vector<double> enb_emissivity_line, int number_sensor, uint32 width_band, vector<double> &surface_temperature_line);
void short_wave_radiation_function(Reader tal_reader, MTL mtl, uint32 width_band, vector<double> &short_wave_radiation_line);
void large_wave_radiation_surface_function(vector<double> eo_emissivity_line, vector<double> surface_temperature_line, uint32 width_band, vector<double> &large_wave_radiation_surface_line);
void large_wave_radiation_atmosphere_function(vector<double> ea_emissivity_line, uint32 width_band, double temperature, vector<double> &large_wave_radiation_atmosphere_line);
void net_radiation_function(vector<double> short_wave_radiation_line, vector<double> large_wave_radiation_surface_line, vector<double> large_wave_radiation_atmosphere_line, vector<double> albedo_line, vector<double> eo_emissivity_line, uint32 width_band, vector<double> &net_radiation_line);
void soil_heat_flux_function(vector<double> ndvi_line, vector<double> surface_temperature_line, vector<double> albedo_line, vector<double> net_radiation_line, uint32 width_band, vector<double> &soil_heat_line);


void d0_fuction(vector<double> pai_line, int width_band, vector<double> &d0_line);
void kb_function(vector<double> ustar_line, vector<double> zom_line, vector<double> pai_line, vector<double> ndvi_line, double ndvi_max, double ndvi_min, int width_band, vector<double> &kb1_line);

void zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, vector<double> d0_line, vector<double> pai_line, int width_band, vector<double> &zom_line);
void zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, int width_band, vector<double> &zom_line);

void ustar_fuction(double u10, vector<double> zom_line, vector<double> d0_line, int width_band, vector<double> &ustar_line);
void ustar_fuction(double u200, vector<double> zom_line, int width_band, vector<double> &ustar_line);

void aerodynamic_resistance_fuction(vector<double> ustar_line, int width_band, vector<double> &aerodynamic_resistance_line);
void aerodynamic_resistance_fuction(vector<double> ustar_line, vector<double> zom_line, vector<double> d0_line, vector<double> kb1_line, int width_band, vector<double> &aerodynamic_resistance_line);

void sensible_heat_function_STEEP(Candidate hot_pixel, Candidate cold_pixel, Station station, uint32 height_band, uint32 width_band, int threads_num, vector<vector<double>> ndvi_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> pai_vector, vector<vector<double>> &sensible_heat_flux_vector);
void sensible_heat_function_default(Candidate hot_pixel, Candidate cold_pixel, Station station, uint32 height_band, uint32 width_band, int threads_num, vector<vector<double>> ndvi_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> &sensible_heat_flux_vector);

void correctionCycleSTEEP(int start_col, int end_col, double a, double b, vector<double> surface_temperature_vector_line,
                        vector<double> d0_vector_line, vector<double> zom_vector_line, vector<double> kb1_vector_line, 
                        vector<double> pai_vector_line, vector<double> ustar_previous_line, vector<double> aerodynamic_resistance_previous_line, 
                        vector<double> &ustar_vector_line, vector<double> &aerodynamic_resistance_vector_line, vector<double> &sensible_heat_flux_vector_line);

void latent_heat_flux_function(vector<double> net_radiation_line, vector<double> soil_heat_flux_line, vector<double> sensible_heat_flux_line, int width_band, vector<double> &latent_heat_flux);
void net_radiation_24h_function(vector<double> albedo_line, double Ra24h, double Rs24h, int width_band, vector<double> &net_radiation_24h_line);
void evapotranspiration_fraction_fuction(vector<double> latent_heat_flux_line, vector<double> net_radiation_line, vector<double> soil_heat_line, int width_band, vector<double> &evapotranspiration_fraction_line);
void sensible_heat_flux_24h_fuction(vector<double> evapotranspiration_fraction_line, vector<double> net_radiation_24h_line, int width_band, vector<double> &sensible_heat_flux_24h_line);
void latent_heat_flux_24h_function(vector<double> evapotranspiration_fraction_line, vector<double> net_radiation_24h_line, int width_band, vector<double> &latent_heat_flux_24h_line);
void evapotranspiration_24h_function(vector<double> latent_heat_flux_24h_line, Station station, int width_band, vector<double> &evapotranspiration_24h_line);
void evapotranspiration_function(vector<double> net_radiation_24h_line, vector<double> evapotranspiration_fraction_line, int width_band, vector<double> &evapotranspiration_line);
