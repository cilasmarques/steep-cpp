#include "landsat.h"

Landsat::Landsat(int method, string bands_paths[], string tal_path, string land_cover_path, int threads_num)
{
  for (int i = 0; i < 8; i++)
    this->bands_paths[i] = bands_paths[i];

  this->tal_path = tal_path;
  this->land_cover_path = land_cover_path;
  this->method = method;
  this->threads_num = threads_num;
};

void Landsat::process_products(MTL mtl, Sensor sensor, Station station)
{
  using namespace std::chrono;
  int64_t general_time, initial_time, final_time, phase2_initial_time, phase2_final_time;
  system_clock::time_point begin, end, phase1_begin, phase2_begin, phase1_end, phase2_end;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

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

  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "Tempo de abertura dos arquivos," << general_time << "," << initial_time << "," << final_time << std::endl;

  // ===== PHASE 1 =====

  phase1_begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  // ==== INITIAL PRODUCTS
  Products products = Products(width_band, height_band);

  for (int line = 0; line < height_band; line++)
  {
    // ==== TAL READER

    system_clock::time_point tal_begin = system_clock::now();
    int64_t tal_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

    tdata_t tal_line_buff = _TIFFmalloc(TIFFScanlineSize(tal));
    unsigned short curr_tal_line_size = TIFFScanlineSize(tal) / width_band;
    Reader tal_reader = Reader(sample_bands, curr_tal_line_size, tal_line_buff);
    TIFFReadScanline(tal, tal_line_buff, line);

    std::cout << "Tempo de leitura do tal," << duration_cast<nanoseconds>(system_clock::now() - tal_begin).count() << "," 
    << tal_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== RADIANCE 

    system_clock::time_point rad_begin = system_clock::now();
    int64_t rad_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.radiance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line);
    std::cout << "Tempo do calculo do rad," << duration_cast<nanoseconds>(system_clock::now() - rad_begin).count() << "," 
    << rad_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== REFLECTANCE

    system_clock::time_point ref_begin = system_clock::now();
    int64_t ref_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.reflectance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line);
    std::cout << "Tempo do calculo do ref," << duration_cast<nanoseconds>(system_clock::now() - ref_begin).count() << ","
    << ref_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== ALBEDO

    system_clock::time_point alb_begin = system_clock::now();
    int64_t alb_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.albedo_function(tal_reader, sensor, width_band, mtl.number_sensor, line);
    std::cout << "Tempo do calculo do alb," << duration_cast<nanoseconds>(system_clock::now() - alb_begin).count() << ","
    << alb_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== NDVI
    system_clock::time_point ndvi_begin = system_clock::now();
    int64_t ndvi_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.ndvi_function(width_band, line);
    std::cout << "Tempo do calculo do ndvi," << duration_cast<nanoseconds>(system_clock::now() - ndvi_begin).count() << ","
    << ndvi_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== LAI
    system_clock::time_point lai_begin = system_clock::now();
    int64_t lai_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.lai_function(width_band, line);
    std::cout << "Tempo do calculo do lai," << duration_cast<nanoseconds>(system_clock::now() - lai_begin).count() << ","
    << lai_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== PAI
    system_clock::time_point pai_begin = system_clock::now();
    int64_t pai_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.pai_function(width_band, line);
    std::cout << "Tempo do calculo do pai," << duration_cast<nanoseconds>(system_clock::now() - pai_begin).count() << ","
    << pai_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== Emissivity indices
    system_clock::time_point emi_begin = system_clock::now();
    int64_t emi_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.enb_emissivity_function(width_band, line);
    products.eo_emissivity_function(width_band, line);
    products.ea_emissivity_function(tal_reader, width_band, line);
    products.surface_temperature_function(mtl.number_sensor, width_band, line);
    std::cout << "Tempo do calculo do TS," << duration_cast<nanoseconds>(system_clock::now() - emi_begin).count() << ","
    << emi_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== Radiation waves
    system_clock::time_point rad_wave_begin = system_clock::now();
    int64_t rad_wave_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.short_wave_radiation_function(tal_reader, mtl, width_band, line);
    products.large_wave_radiation_surface_function(width_band, line);
    products.large_wave_radiation_atmosphere_function(width_band, station.temperature_image, line);
    std::cout << "Tempo do calculo do rad_wave," << duration_cast<nanoseconds>(system_clock::now() - rad_wave_begin).count() << ","
    << rad_wave_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== Rn
    system_clock::time_point rn_begin = system_clock::now();
    int64_t rn_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.net_radiation_function(width_band, line);
    std::cout << "Tempo do calculo do rn," << duration_cast<nanoseconds>(system_clock::now() - rn_begin).count() << ","
    << rn_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    // ==== G
    system_clock::time_point g_begin = system_clock::now();
    int64_t g_init = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    products.soil_heat_flux_function(width_band, line);
    std::cout << "Tempo do calculo do g," << duration_cast<nanoseconds>(system_clock::now() - g_begin).count() << ","
    << g_init << "," << duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count() << std::endl;

    _TIFFfree(tal_line_buff);
  }
  phase1_end = system_clock::now();
  general_time = duration_cast<nanoseconds>(phase1_end - phase1_begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P1 - TOTAL," << general_time << "," << initial_time << "," << final_time << std::endl;

  TIFFClose(tal);

  // ===== PHASE 2 =====

  phase2_begin = system_clock::now();
  phase2_initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  // ==== PIXEL SELECTION

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  Candidate hot_pixel, cold_pixel;
  if (this->method == 0)
  { // STEEP
    pair<Candidate, Candidate> pixels = getColdHotPixelsSTEPP(products.ndvi_vector, products.surface_temperature_vector, products.albedo_vector, products.net_radiation_vector, products.soil_heat_vector, height_band, width_band);
    hot_pixel = pixels.first, cold_pixel = pixels.second;    
  }
  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - PIXEL SELECTION," << general_time << "," << initial_time << "," << final_time << std::endl;

  // ==== RAH CYCLE - COMPUTE H

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  if (this->method == 0)
  { // STEEP
    products.sensible_heat_function_STEEP(hot_pixel, cold_pixel, station, height_band, width_band, this->threads_num);
  }
  end = system_clock::now();
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - H," << general_time << "," << initial_time << "," << final_time << std::endl;

  // ==== FINAL PRODUCTS

  double dr = (1 / mtl.distance_earth_sun) * (1 / mtl.distance_earth_sun);
  double sigma = 0.409 * sin(((2 * PI / 365) * mtl.julian_day) - 1.39);
  double phi = (PI / 180) * station.latitude;
  double omegas = acos(-tan(phi) * tan(sigma));
  double Ra24h = (((24 * 60 / PI) * GSC * dr) * (omegas * sin(phi) * sin(sigma) + cos(phi) * cos(sigma) * sin(omegas))) * (1000000 / 86400.0);
  double Rs24h = station.INTERNALIZATION_FACTOR * sqrt(station.v7_max - station.v7_min) * Ra24h;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
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
  general_time = duration_cast<nanoseconds>(end - begin).count();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - FINAL PRODUCTS," << general_time << "," << initial_time << "," << final_time << std::endl;

  phase2_end = system_clock::now();
  general_time = duration_cast<nanoseconds>(phase2_end - phase2_begin).count();
  phase2_final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  std::cout << "P2 - TOTAL," << general_time << "," << phase2_initial_time << "," << phase2_final_time << std::endl;

  // =====  END + OUTPUTS =====

  string outPath = "./output/extra.txt";
  std::ofstream outputThreads(outPath);
  std::streambuf* coutThreads = std::cout.rdbuf();
  std::cout.rdbuf(outputThreads.rdbuf());
  std::cout << "threads_num: " << threads_num << std::endl;
  std::cout << "width: " << width_band << std::endl;
  std::cout << "height: " << height_band << std::endl;

  // std::ofstream outputProds("./output/products.txt");
  // std::streambuf* coutProds = std::cout.rdbuf();
  // std::cout.rdbuf(outputProds.rdbuf());

  // std::cout << " ==== albedo" << std::endl;
  // printVector2x2(products.albedo_vector);

  // std::cout << " ==== ndvi" << std::endl;
  // printVector2x2(products.ndvi_vector);

  // std::cout << " ==== net_radiation" << std::endl;
  // printVector2x2(products.net_radiation_vector);

  // std::cout << " ==== soil_heat" << std::endl;
  // printVector2x2(products.soil_heat_vector);

  // std::cout << " ==== sensible_heat_flux" << std::endl;
  // printVector2x2(products.sensible_heat_flux_vector);

  // std::cout << " ==== latent_heat_flux" << std::endl;
  // printVector2x2(products.latent_heat_flux_vector);

  // std::cout << " ==== net_radiation_24h" << std::endl;
  // printVector2x2(products.net_radiation_24h_vector);

  // std::cout << " ==== evapotranspiration_fraction" << std::endl;
  // printVector2x2(products.evapotranspiration_fraction_vector);

  // std::cout << " ==== sensible_heat_flux_24h" << std::endl;
  // printVector2x2(products.sensible_heat_flux_24h_vector);

  // std::cout << " ==== latent_heat_flux_24h" << std::endl;
  // printVector2x2(products.latent_heat_flux_24h_vector);

  // std::cout << " ==== evapotranspiration_24h" << std::endl;
  // printVector2x2(products.evapotranspiration_24h_vector);

  // std::cout << " ==== evapotranspiration" << std::endl;
  // printVector2x2(products.evapotranspiration_vector);
};
