#include "landsat.h"

Landsat::Landsat(int method, string bands_paths[], string tal_path, string land_cover_path)
{
  for (int i = 0; i < 8; i++)
    this->bands_paths[i] = bands_paths[i];

  this->tal_path = tal_path;
  this->land_cover_path = land_cover_path;
  this->method = method;
};

void Landsat::process_products(MTL mtl, Sensor sensor, Station station)
{
  Reader TIFFs_reader = Reader();

  TIFF *bands_resampled[8];
  for (int i = 1; i < 8; i++)
  {
    string path_tiff_base = this->bands_paths[i];
    bands_resampled[i] = TIFFOpen(path_tiff_base.c_str(), "rm");
    TIFFs_reader.check_open_tiff(bands_resampled[i]);
  }

  TIFF *tal = TIFFOpen(this->tal_path.c_str(), "rm");
  TIFFs_reader.check_open_tiff(tal);

  uint16 sample_bands;
  uint32 height_band, width_band;
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGELENGTH, &height_band);
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGEWIDTH, &width_band);
  TIFFGetField(bands_resampled[1], TIFFTAG_SAMPLEFORMAT, &sample_bands);

  // Declare auxiliaries arrays
  vector<vector<vector<double>>> radiance_vector(height_band, vector<vector<double>>(width_band, vector<double>(8)));
  vector<vector<vector<double>>> reflectance_vector(height_band, vector<vector<double>>(width_band, vector<double>(8)));

  vector<vector<double>> albedo_vector(height_band, vector<double>(width_band));
  vector<vector<double>> ndvi_vector(height_band, vector<double>(width_band));
  vector<vector<double>> soil_heat_vector(height_band, vector<double>(width_band));
  vector<vector<double>> surface_temperature_vector(height_band, vector<double>(width_band));
  vector<vector<double>> net_radiation_vector(height_band, vector<double>(width_band));
  vector<vector<double>> lai_vector(height_band, vector<double>(width_band));
  vector<vector<double>> evi_vector(height_band, vector<double>(width_band));
  vector<vector<double>> pai_vector(height_band, vector<double>(width_band));
  vector<vector<double>> enb_emissivity_vector(height_band, vector<double>(width_band));
  vector<vector<double>> eo_emissivity_vector(height_band, vector<double>(width_band));
  vector<vector<double>> ea_emissivity_vector(height_band, vector<double>(width_band));
  vector<vector<double>> short_wave_radiation_vector(height_band, vector<double>(width_band));
  vector<vector<double>> large_wave_radiation_surface_vector(height_band, vector<double>(width_band));
  vector<vector<double>> large_wave_radiation_atmosphere_vector(height_band, vector<double>(width_band));

  for (int line = 0; line < height_band; line++)
  {
    tdata_t tal_line_buff = _TIFFmalloc(TIFFScanlineSize(tal));
    unsigned short curr_tal_line_size = TIFFScanlineSize(tal) / width_band;
    Reader tal_reader = Reader(sample_bands, curr_tal_line_size, tal_line_buff);
    TIFFReadScanline(tal, tal_line_buff, line);

    radiance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line, radiance_vector[line]);
    reflectance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line, reflectance_vector[line]);
    albedo_function(tal_reader, sensor, width_band, mtl.number_sensor, reflectance_vector[line], albedo_vector[line]);

    // Vegetation indices
    ndvi_function(reflectance_vector[line], width_band, ndvi_vector[line]);
    pai_function(reflectance_vector[line], width_band, pai_vector[line]);
    lai_function(reflectance_vector[line], width_band, lai_vector[line]);
    evi_function(reflectance_vector[line], width_band, evi_vector[line]);

    // Emissivity indices
    ea_emissivity_function(tal_reader, width_band, ea_emissivity_vector[line]);
    enb_emissivity_function(lai_vector[line], ndvi_vector[line], width_band, enb_emissivity_vector[line]);
    eo_emissivity_function(lai_vector[line], ndvi_vector[line], width_band, eo_emissivity_vector[line]);
    surface_temperature_function(radiance_vector[line], enb_emissivity_vector[line], mtl.number_sensor, width_band, surface_temperature_vector[line]);

    // Radiation waves
    short_wave_radiation_function(tal_reader, mtl, width_band, short_wave_radiation_vector[line]);
    large_wave_radiation_surface_function(eo_emissivity_vector[line], surface_temperature_vector[line], width_band, large_wave_radiation_surface_vector[line]);
    large_wave_radiation_atmosphere_function(ea_emissivity_vector[line], width_band, station.temperature_image, large_wave_radiation_atmosphere_vector[line]);

    // Main products
    net_radiation_function(short_wave_radiation_vector[line], large_wave_radiation_surface_vector[line], large_wave_radiation_atmosphere_vector[line], albedo_vector[line], eo_emissivity_vector[line], width_band, net_radiation_vector[line]);
    soil_heat_flux_function(ndvi_vector[line], surface_temperature_vector[line], albedo_vector[line], net_radiation_vector[line], width_band, soil_heat_vector[line]);

    _TIFFfree(tal_line_buff);
  }

  TIFFClose(tal);

  Candidate hot_pixel, cold_pixel;
  if (this->method == 0)
  { // STEEP
    hot_pixel = getHotPixelSTEPP(ndvi_vector, surface_temperature_vector, albedo_vector, net_radiation_vector, soil_heat_vector, height_band, width_band);
    cold_pixel = getColdPixelSTEPP(ndvi_vector, surface_temperature_vector, albedo_vector, net_radiation_vector, soil_heat_vector, height_band, width_band);
  }
  else if (this->method == 1)
  { // ASEBAL
    hot_pixel = getHotPixelASEBAL(ndvi_vector, surface_temperature_vector, albedo_vector, net_radiation_vector, soil_heat_vector, height_band, width_band);
    cold_pixel = getColdPixelASEBAL(ndvi_vector, surface_temperature_vector, albedo_vector, net_radiation_vector, soil_heat_vector, height_band, width_band);
  }
  // else if (this->method == 2)
  // { // ESA SEBAL
  //   TIFF *land_cover = TIFFOpen(this->land_cover_path.c_str(), "r");
  //   pair<Candidate, Candidate> pixels = getColdHotPixelsESA(&ndvi, &surface_temperature, &albedo, &net_radiation, &soil_heat, &land_cover, height_band, width_band, this->output_path);
  //   hot_pixel = pixels.first, cold_pixel = pixels.second;
  // }

  vector<vector<double>> sensible_heat_flux_vector(height_band, vector<double>(width_band));
  if (this->method == 0)
  { // STEEP
    sensible_heat_function_STEEP(hot_pixel, cold_pixel, station, height_band, width_band, ndvi_vector, net_radiation_vector, soil_heat_vector, surface_temperature_vector, pai_vector, sensible_heat_flux_vector);
  }
  else
  {
    sensible_heat_function_default(hot_pixel, cold_pixel, station, height_band, width_band, ndvi_vector, net_radiation_vector, soil_heat_vector, surface_temperature_vector, sensible_heat_flux_vector);
  }

  vector<vector<double>> latent_heat_flux_vector(height_band, vector<double>(width_band));
  vector<vector<double>> net_radiation_24h_vector(height_band, vector<double>(width_band));
  vector<vector<double>> evapotranspiration_fraction_vector(height_band, vector<double>(width_band));
  vector<vector<double>> sensible_heat_flux_24h_vector(height_band, vector<double>(width_band));
  vector<vector<double>> latent_heat_flux_24h_vector(height_band, vector<double>(width_band));
  vector<vector<double>> evapotranspiration_24h_vector(height_band, vector<double>(width_band));
  vector<vector<double>> evapotranspiration_vector(height_band, vector<double>(width_band));

  // Upscalling temporal
  double dr = (1 / mtl.distance_earth_sun) * (1 / mtl.distance_earth_sun);
  double sigma = 0.409 * sin(((2 * PI / 365) * mtl.julian_day) - 1.39);
  double phi = (PI / 180) * station.latitude;
  double omegas = acos(-tan(phi) * tan(sigma));
  double Ra24h = (((24 * 60 / PI) * GSC * dr) * (omegas * sin(phi) * sin(sigma) + cos(phi) * cos(sigma) * sin(omegas))) * (1000000 / 86400.0);

  // Short wave radiation incident in 24 hours (Rs24h)
  double Rs24h = station.INTERNALIZATION_FACTOR * sqrt(station.v7_max - station.v7_min) * Ra24h;

  for (int line = 0; line < height_band; line++)
  {
    latent_heat_flux_function(net_radiation_vector[line], soil_heat_vector[line], sensible_heat_flux_vector[line], width_band, latent_heat_flux_vector[line]);
    net_radiation_24h_function(albedo_vector[line], Ra24h, Rs24h, width_band, net_radiation_24h_vector[line]);
    evapotranspiration_fraction_fuction(latent_heat_flux_vector[line], net_radiation_vector[line], soil_heat_vector[line], width_band, evapotranspiration_fraction_vector[line]);
    sensible_heat_flux_24h_fuction(evapotranspiration_fraction_vector[line], net_radiation_24h_vector[line], width_band, sensible_heat_flux_24h_vector[line]);
    latent_heat_flux_24h_function(evapotranspiration_fraction_vector[line], net_radiation_24h_vector[line], width_band, latent_heat_flux_24h_vector[line]);
    evapotranspiration_24h_function(latent_heat_flux_24h_vector[line], station, width_band, evapotranspiration_24h_vector[line]);
    evapotranspiration_function(net_radiation_24h_vector[line], evapotranspiration_fraction_vector[line], width_band, evapotranspiration_vector[line]);
  }

  std::cout << " ==== albedo" << std::endl;
  printVector2x2(albedo_vector);

  std::cout << " ==== ndvi" << std::endl;
  printVector2x2(ndvi_vector);

  std::cout << " ==== net_radiation" << std::endl;
  printVector2x2(net_radiation_vector);

  std::cout << " ==== soil_heat" << std::endl;
  printVector2x2(soil_heat_vector);

  std::cout << " ==== sensible_heat_flux" << std::endl;
  printVector2x2(sensible_heat_flux_vector);

  std::cout << " ==== latent_heat_flux" << std::endl;
  printVector2x2(latent_heat_flux_vector);

  std::cout << " ==== net_radiation_24h" << std::endl;
  printVector2x2(net_radiation_24h_vector);

  std::cout << " ==== evapotranspiration_fraction" << std::endl;
  printVector2x2(evapotranspiration_fraction_vector);

  std::cout << " ==== sensible_heat_flux_24h" << std::endl;
  printVector2x2(sensible_heat_flux_24h_vector);

  std::cout << " ==== latent_heat_flux_24h" << std::endl;
  printVector2x2(latent_heat_flux_24h_vector);

  std::cout << " ==== evapotranspiration_24h" << std::endl;
  printVector2x2(evapotranspiration_24h_vector);

  std::cout << " ==== evapotranspiration" << std::endl;
  printVector2x2(evapotranspiration_vector);
};
