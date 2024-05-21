#include "landsat.h"

Landsat::Landsat(string bands_paths[], string land_cover_path, MTL mtl)
{
  Reader TIFFs_reader = Reader();

  for (int i = 1; i <= 8; i++)
  {
    string path_tiff_base = bands_paths[i];
    this->bands_resampled[i] = TIFFOpen(path_tiff_base.c_str(), "rm");
    TIFFs_reader.check_open_tiff(this->bands_resampled[i]);
  }

  uint16_t sample_format;
  uint32_t height, width;
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGELENGTH, &height);
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(bands_resampled[1], TIFFTAG_SAMPLEFORMAT, &sample_format);

  this->width_band = width;
  this->height_band = height;
  this->sample_bands = sample_format;

  this->mtl = mtl;
  this->products = Products(this->width_band, this->height_band);
};

string Landsat::select_endmembers(int method)
{
  system_clock::time_point begin, end;
  int64_t general_time, initial_time, final_time;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  if (method == 0)
  { // STEEP
	  pair<Candidate, Candidate> pixels = getEndmembersSTEPP(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
		hot_pixel = pixels.first; cold_pixel = pixels.second;
  }
  else if (method == 1)
  { // ASEBAL
	  pair<Candidate, Candidate> pixels = getEndmembersASEBAL(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
		hot_pixel = pixels.first; cold_pixel = pixels.second;
  }

  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  return "P2 - PIXEL SELECTION," + std::to_string(general_time) + "," + std::to_string(initial_time) + "," + std::to_string(final_time) + "\n";
}

string Landsat::converge_rah_cycle(Station station, int method, int threads_num, int blocks_num)
{
  string result = "";
  system_clock::time_point begin, end;
  int64_t general_time, initial_time, final_time;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  double ustar_station = (VON_KARMAN * station.v6) / (log(station.WIND_SPEED / station.SURFACE_ROUGHNESS));
  double u10 = (ustar_station / VON_KARMAN) * log(10 / station.SURFACE_ROUGHNESS);
  double ndvi_min = 1.0;
  double ndvi_max = -1.0;

  for (int line = 0; line < this->height_band; line++)
  {
    vector<double> ndvi_line = products.ndvi_vector[line];
    for (int col = 0; col < this->width_band; col++)
    {
      if (ndvi_line[col] < ndvi_min)
        ndvi_min = ndvi_line[col];
      if (ndvi_line[col] > ndvi_max)
        ndvi_max = ndvi_line[col];
    }
  }

  for (int line = 0; line < height_band; line++)
  {
    products.d0_fuction(line);
    products.zom_fuction(station.A_ZOM, station.B_ZOM, line);
    products.ustar_fuction(u10, line);
    products.kb_function(ndvi_max, ndvi_min, line);
    products.aerodynamic_resistance_fuction(line);
  }

  result += products.rah_correction_function_blocks(ndvi_min, ndvi_max, hot_pixel, cold_pixel);

  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  result += "P2 - RAH CYCLE," + std::to_string(general_time) + "," + std::to_string(initial_time) + "," + std::to_string(final_time) + "\n";
  return result;
};


string Landsat::compute_Rn_G(Sensor sensor, Station station)
{
  system_clock::time_point begin, end;
  int64_t general_time, initial_time, final_time;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  TIFF *tal = this->bands_resampled[8];
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
  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  return "P1 - Rn_G," + std::to_string(general_time) + "," + std::to_string(initial_time) + "," + std::to_string(final_time) + "\n";
}

string Landsat::compute_H_ET(Station station)
{
  system_clock::time_point begin, end;
  int64_t general_time, initial_time, final_time;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  double dr = (1 / mtl.distance_earth_sun) * (1 / mtl.distance_earth_sun);
  double sigma = 0.409 * sin(((2 * PI / 365) * mtl.julian_day) - 1.39);
  double phi = (PI / 180) * station.latitude;
  double omegas = acos(-tan(phi) * tan(sigma));
  double Ra24h = (((24 * 60 / PI) * GSC * dr) * (omegas * sin(phi) * sin(sigma) + cos(phi) * cos(sigma) * sin(omegas))) * (1000000 / 86400.0);
  double Rs24h = station.INTERNALIZATION_FACTOR * sqrt(station.v7_max - station.v7_min) * Ra24h;

  double dt_pq_terra = products.H_pq_terra * products.rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);
  double dt_pf_terra = products.H_pf_terra * products.rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);

  double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
  double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

  for (int line = 0; line < height_band; line++)
  {
    products.sensible_heat_flux_function(a, b, line);
    products.latent_heat_flux_function(width_band, line);
    products.net_radiation_24h_function(Ra24h, Rs24h, width_band, line);
    products.evapotranspiration_fraction_fuction(width_band, line);
    products.sensible_heat_flux_24h_fuction(width_band, line);
    products.latent_heat_flux_24h_function(width_band, line);
    products.evapotranspiration_24h_function(station, width_band, line);
    products.evapotranspiration_function(width_band, line);
  }
  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  return "P2 - FINAL PRODUCTS," + std::to_string(general_time) + "," + std::to_string(initial_time) + "," + std::to_string(final_time) + "\n";
};

void Landsat::save_products(string output_path)
{
  system_clock::time_point begin, end;
  int64_t general_time, initial_time, final_time;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  std::ofstream outputProds(output_path);
  std::streambuf *coutProds = std::cout.rdbuf();
  std::cout.rdbuf(outputProds.rdbuf());

  std::cout << " ==== albedo" << std::endl;
  printVector2x2(products.albedo_vector);

  std::cout << " ==== ndvi" << std::endl;
  printVector2x2(products.ndvi_vector);

  std::cout << " ==== net_radiation" << std::endl;
  printVector2x2(products.net_radiation_vector);

  std::cout << " ==== soil_heat" << std::endl;
  printVector2x2(products.soil_heat_vector);

  std::cout << " ==== sensible_heat_flux" << std::endl;
  printVector2x2(products.sensible_heat_flux_vector);

  std::cout << " ==== latent_heat_flux" << std::endl;
  printVector2x2(products.latent_heat_flux_vector);

  std::cout << " ==== net_radiation_24h" << std::endl;
  printVector2x2(products.net_radiation_24h_vector);

  std::cout << " ==== evapotranspiration_fraction" << std::endl;
  printVector2x2(products.evapotranspiration_fraction_vector);

  std::cout << " ==== sensible_heat_flux_24h" << std::endl;
  printVector2x2(products.sensible_heat_flux_24h_vector);

  std::cout << " ==== latent_heat_flux_24h" << std::endl;
  printVector2x2(products.latent_heat_flux_24h_vector);

  std::cout << " ==== evapotranspiration_24h" << std::endl;
  printVector2x2(products.evapotranspiration_24h_vector);

  std::cout << " ==== evapotranspiration" << std::endl;
  printVector2x2(products.evapotranspiration_vector);

  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P3 - WRITE PRODUCTS," << general_time << "," << initial_time << "," << final_time << std::endl;
};

void Landsat::close()
{
  for (int i = 1; i <= 8; i++)
  {
    TIFFClose(this->bands_resampled[i]);
  }
};