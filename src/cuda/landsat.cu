#include "landsat.h"

Landsat::Landsat(int method, string bands_paths[], string tal_path, string land_cover_path, int threads_num, int blocks_num)
{
  for (int i = 0; i < 8; i++)
    this->bands_paths[i] = bands_paths[i];

  this->tal_path = tal_path;
  this->land_cover_path = land_cover_path;
  this->method = method;
  this->threads_num = threads_num;
  this->blocks_num = blocks_num;
};

void Landsat::process_products(MTL mtl, Sensor sensor, Station station)
{
  using namespace std::chrono;
  int64_t general_time, initial_time, final_time, phase2_initial_time, phase2_final_time;
  system_clock::time_point begin, end, phase1_begin, phase2_begin, phase1_end, phase2_end;

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

  uint16_t sample_bands;
  uint32_t height_band, width_band;
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGELENGTH, &height_band);
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGEWIDTH, &width_band);
  TIFFGetField(bands_resampled[1], TIFFTAG_SAMPLEFORMAT, &sample_bands);

  // ===== PHASE 1 =====

  phase1_begin = system_clock::now();
  initial_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

  // ==== INITIAL PRODUCTS
  Products products = Products(width_band, height_band);

  for (int line = 0; line < height_band; line++)
  {
    tdata_t tal_line_buff = _TIFFmalloc(TIFFScanlineSize(tal));
    unsigned short curr_tal_line_size = TIFFScanlineSize(tal) / width_band;
    Reader tal_reader = Reader(sample_bands, curr_tal_line_size, tal_line_buff);
    TIFFReadScanline(tal, tal_line_buff, line);

    products.radiance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line);
    products.reflectance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line);
    products.albedo_function(tal_reader, sensor, width_band, mtl.number_sensor, line);

    // Vegetation indices
    products.ndvi_function(width_band, line);
    products.pai_function(width_band, line);
    products.lai_function(width_band, line);
    products.evi_function(width_band, line);

    // Emissivity indices
    products.enb_emissivity_function(width_band, line);
    products.eo_emissivity_function(width_band, line);
    products.ea_emissivity_function(tal_reader, width_band, line);
    products.surface_temperature_function(mtl.number_sensor, width_band, line);

    // Radiation waves
    products.short_wave_radiation_function(tal_reader, mtl, width_band, line);
    products.large_wave_radiation_surface_function(width_band, line);
    products.large_wave_radiation_atmosphere_function(width_band, station.temperature_image, line);

    // Main products
    products.net_radiation_function(width_band, line);
    products.soil_heat_flux_function(width_band, line);

    _TIFFfree(tal_line_buff);
  }
  phase1_end = system_clock::now();
  general_time = duration_cast<milliseconds>(phase1_end - phase1_begin).count();
  final_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P1 - TOTAL," << general_time << "," << initial_time << "," << final_time << std::endl;

  TIFFClose(tal);

  // ===== PHASE 2 =====

  phase2_begin = system_clock::now();
  phase2_initial_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

  // ==== PIXEL SELECTION

  begin = system_clock::now();
  initial_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  Candidate hot_pixel, cold_pixel;
  if (this->method == 0)
  { // STEEP
    hot_pixel = getHotPixelSTEPP(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
    cold_pixel = getColdPixelSTEPP(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
  }
  else if (this->method == 1)
  { // ASEBAL
    hot_pixel = getHotPixelASEBAL(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
    cold_pixel = getColdPixelASEBAL(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
  }
  else if (this->method == 2)
  { // ESA SEBAL
    TIFF *land_cover_tiff = TIFFOpen(this->land_cover_path.c_str(), "r");
    pair<Candidate, Candidate> pixels = getColdHotPixelsESA(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, land_cover_tiff, height_band, width_band);
    hot_pixel = pixels.first, cold_pixel = pixels.second;
  }
  end = system_clock::now();
  general_time = duration_cast<milliseconds>(end - begin).count();
  final_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - PIXEL SELECTION," << general_time << "," << initial_time << "," << final_time << std::endl;

  // ==== RAH CYCLE - COMPUTE H

  begin = system_clock::now();
  initial_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  if (this->method == 0)
  { // STEEP
    products.sensible_heat_function_STEEP(hot_pixel, cold_pixel, station, height_band, width_band, this->threads_num, this->blocks_num);
  }
  // else
  // { // ASEBAL & ESASEB
  //   sensible_heat_function_default(hot_pixel, cold_pixel, station, height_band, width_band, this->threads_num, this->blocks_num);
  // }
  end = system_clock::now();
  general_time = duration_cast<milliseconds>(end - begin).count();
  final_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - H," << general_time << "," << initial_time << "," << final_time << std::endl;

  // ==== FINAL PRODUCTS

  double dr = (1 / mtl.distance_earth_sun) * (1 / mtl.distance_earth_sun);
  double sigma = 0.409 * sin(((2 * PI / 365) * mtl.julian_day) - 1.39);
  double phi = (PI / 180) * station.latitude;
  double omegas = acos(-tan(phi) * tan(sigma));
  double Ra24h = (((24 * 60 / PI) * GSC * dr) * (omegas * sin(phi) * sin(sigma) + cos(phi) * cos(sigma) * sin(omegas))) * (1000000 / 86400.0);
  double Rs24h = station.INTERNALIZATION_FACTOR * sqrt(station.v7_max - station.v7_min) * Ra24h;

  begin = system_clock::now();
  initial_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  for (int line = 0; line < height_band; line++)
  {
    products.latent_heat_flux_function(width_band, line);
    products.net_radiation_24h_function(Ra24h, Rs24h, width_band, line);
    products.evapotranspiration_fraction_fuction(width_band, line);
    products.sensible_heat_flux_24h_fuction(width_band, line);
    products.latent_heat_flux_24h_function(width_band, line);
    products.evapotranspiration_24h_function(station, width_band, line);
    products.evapotranspiration_function(width_band, line);
  }
  end = system_clock::now();
  general_time = duration_cast<milliseconds>(end - begin).count();
  final_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - FINAL PRODUCTS," << general_time << "," << initial_time << "," << final_time << std::endl;

  phase2_end = system_clock::now();
  general_time = duration_cast<milliseconds>(phase2_end - phase2_begin).count();
  phase2_final_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - TOTAL," << general_time << "," << phase2_initial_time << "," << phase2_final_time << std::endl;
};
