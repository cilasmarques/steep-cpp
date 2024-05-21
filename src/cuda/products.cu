#include "products.h"
#include "kernels.cuh"
#include "cuda_utils.h"

Products::Products() {}

Products::Products(uint32_t width_band, uint32_t height_band)
{
  this->width_band = width_band;
  this->height_band = height_band;
  this->nBytes_band = height_band * width_band * sizeof(double);

  this->radiance_vector = vector<vector<vector<double>>>(height_band, vector<vector<double>>(width_band, vector<double>(8)));
  this->reflectance_vector = vector<vector<vector<double>>>(height_band, vector<vector<double>>(width_band, vector<double>(8)));
  this->albedo_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->ndvi_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->soil_heat_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->surface_temperature_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->net_radiation_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->lai_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->evi_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->pai_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->enb_emissivity_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->eo_emissivity_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->ea_emissivity_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->short_wave_radiation_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->large_wave_radiation_surface_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->large_wave_radiation_atmosphere_vector = vector<vector<double>>(height_band, vector<double>(width_band));

  // ==== GPU
  this->surface_temperature_pointer = (double *)malloc(nBytes_band);
  this->d0_pointer = (double *)malloc(nBytes_band);
  this->zom_pointer = (double *)malloc(nBytes_band);
  this->ustar_pointer = (double *)malloc(nBytes_band);
  this->kb1_pointer = (double *)malloc(nBytes_band);
  this->aerodynamic_resistance_pointer = (double *)malloc(nBytes_band);
  this->sensible_heat_flux_pointer = (double *)malloc(nBytes_band);

  this->sensible_heat_flux_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->latent_heat_flux_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->net_radiation_24h_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->evapotranspiration_fraction_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->sensible_heat_flux_24h_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->latent_heat_flux_24h_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->evapotranspiration_24h_vector = vector<vector<double>>(height_band, vector<double>(width_band));
  this->evapotranspiration_vector = vector<vector<double>>(height_band, vector<double>(width_band));
};

void Products::radiance_function(TIFF **bands_resampled, uint32_t width_band, uint16_t sample_bands, MTL mtl, Sensor sensor, int line)
{
  int inital_index = (mtl.number_sensor == 8) ? 7 : 1;

  for (int i_band = inital_index; i_band < 8; i_band++)
  {
    TIFF *curr_band = bands_resampled[i_band];
    tdata_t band_line_buff = _TIFFmalloc(TIFFScanlineSize(curr_band));
    unsigned short curr_band_line_size = TIFFScanlineSize(curr_band) / width_band;
    Reader pixel_reader = Reader(sample_bands, curr_band_line_size, band_line_buff);

    TIFFReadScanline(curr_band, band_line_buff, line);

    for (int col = 0; col < width_band; col++)
    {
      double radiance_pixel = NaN;
      double band_pixel = pixel_reader.read_tiff_pixel(col);

      if (band_pixel > 0)
      {
        // if (mtl.number_sensor == 8)
        //   radiance_pixel = band_pixel * mtl.rad_mult_10 + mtl.rad_add_10;
        // else
        radiance_pixel = band_pixel * sensor.parameters[i_band][sensor.GRESCALE] + sensor.parameters[i_band][sensor.BRESCALE];
      }

      this->radiance_vector[line][col][i_band] = radiance_pixel;
    }

    _TIFFfree(band_line_buff);
  }
}

void Products::reflectance_function(TIFF **bands_resampled, uint32_t width_band, uint16_t sample_bands, MTL mtl, Sensor sensor, int line)
{
  for (int i_band = 1; i_band < 8; i_band++)
  {
    TIFF *curr_band = bands_resampled[i_band];
    tdata_t band_line_buff = _TIFFmalloc(TIFFScanlineSize(curr_band));
    unsigned short curr_band_line_size = TIFFScanlineSize(curr_band) / width_band;
    Reader pixel_reader = Reader(sample_bands, curr_band_line_size, band_line_buff);

    TIFFReadScanline(curr_band, band_line_buff, line);

    for (int col = 0; col < width_band; col++)
    {
      double reflectance_pixel = NaN;
      double band_pixel = pixel_reader.read_tiff_pixel(col);

      if (band_pixel > 0)
      {
        double radiance_pixel = band_pixel * sensor.parameters[i_band][sensor.GRESCALE] + sensor.parameters[i_band][sensor.BRESCALE];

        if (mtl.number_sensor == 8)
          reflectance_pixel = radiance_pixel / sin(mtl.sun_elevation * PI / 180);
        else
          reflectance_pixel = (PI * radiance_pixel * mtl.distance_earth_sun * mtl.distance_earth_sun) / (sensor.parameters[i_band][sensor.ESUN] * sin(mtl.sun_elevation * PI / 180));
      }

      this->reflectance_vector[line][col][i_band] = reflectance_pixel;
    }

    _TIFFfree(band_line_buff);
  }
}

void Products::albedo_function(Reader tal_reader, Sensor sensor, uint32_t width_band, int number_sensor, int line)
{
  int final_tif_calc = number_sensor == 8 ? 6 : 7;

  for (int col = 0; col < width_band; col++)
  {
    double alb = this->reflectance_vector[line][col][1] * sensor.parameters[1][sensor.WB] +
                 this->reflectance_vector[line][col][2] * sensor.parameters[2][sensor.WB] +
                 this->reflectance_vector[line][col][3] * sensor.parameters[3][sensor.WB] +
                 this->reflectance_vector[line][col][4] * sensor.parameters[4][sensor.WB] +
                 this->reflectance_vector[line][col][5] * sensor.parameters[5][sensor.WB] +
                 this->reflectance_vector[line][col][final_tif_calc] * sensor.parameters[final_tif_calc][sensor.WB];

    alb = (alb - 0.03) / (tal_reader.read_tiff_pixel(col) * tal_reader.read_tiff_pixel(col));

    this->albedo_vector[line][col] = alb;
  }
}

void Products::ndvi_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    this->ndvi_vector[line][col] =
        (this->reflectance_vector[line][col][4] - this->reflectance_vector[line][col][3]) /
        (this->reflectance_vector[line][col][4] + this->reflectance_vector[line][col][3]);
  }
};

void Products::pai_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    double pai_value = 10.1 * (this->reflectance_vector[line][col][4] - sqrt(this->reflectance_vector[line][col][3])) + 3.1;

    if (pai_value < 0)
      pai_value = 0;

    this->pai_vector[line][col] = pai_value;
  }
};

void Products::lai_function(uint32_t width_band, int line)
{
  double savi_line[width_band];
  double L = 0.05;

  for (int col = 0; col < width_band; col++)
  {
    savi_line[col] = ((1 + L) * (this->reflectance_vector[line][col][4] - this->reflectance_vector[line][col][3])) /
                     (L + (this->reflectance_vector[line][col][4] + this->reflectance_vector[line][col][3]));

    if (!isnan(savi_line[col]) && definitelyGreaterThan(savi_line[col], 0.687))
      this->lai_vector[line][col] = 6;
    else if (!isnan(savi_line[col]) && definitelyLessThan(savi_line[col], 0.1))
      this->lai_vector[line][col] = 0;
    else
      this->lai_vector[line][col] = -log((0.69 - savi_line[col]) / 0.59) / 0.91;
  }
};

void Products::evi_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    this->evi_vector[line][col] = 2.5 * ((this->reflectance_vector[line][col][4] - this->reflectance_vector[line][col][3]) /
                                         (this->reflectance_vector[line][col][4] + 6 * this->reflectance_vector[line][col][3] - 7.5 * this->reflectance_vector[line][col][1] + 1));
  }
};

void Products::enb_emissivity_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    if (definitelyLessThan(this->ndvi_vector[line][col], 0) || definitelyGreaterThan(this->lai_vector[line][col], 2.99))
      this->enb_emissivity_vector[line][col] = 0.98;
    else
      this->enb_emissivity_vector[line][col] = 0.97 + 0.0033 * this->lai_vector[line][col];
  }
};

void Products::eo_emissivity_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    if (definitelyLessThan(this->ndvi_vector[line][col], 0) || definitelyGreaterThan(this->lai_vector[line][col], 2.99))
      this->eo_emissivity_vector[line][col] = 0.98;
    else
      this->eo_emissivity_vector[line][col] = 0.95 + 0.01 * this->lai_vector[line][col];
  }
};

void Products::ea_emissivity_function(Reader tal_reader, uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
    this->ea_emissivity_vector[line][col] = 0.85 * pow((-1 * log(tal_reader.read_tiff_pixel(col))), 0.09);
};

void Products::surface_temperature_function(int number_sensor, uint32_t width_band, int line)
{
  double k1, k2;

  switch (number_sensor)
  {
  case 5:
    k1 = 607.76;
    k2 = 1282.71;

    break;
  case 7:
    k1 = 666.09;
    k2 = 1260.56;

    break;
  case 8:
    k1 = 774.8853;
    k2 = 1321.0789;

    break;
  default:
    cerr << "Sensor problem!";
    exit(6);
  }

  int radiance_number = (number_sensor == 5) ? 6 : 7;

  for (int col = 0; col < width_band; col++)
  {
    double surface_temperature_vector = k2 / (log((this->enb_emissivity_vector[line][col] * k1 / this->radiance_vector[line][col][radiance_number]) + 1));
    this->surface_temperature_pointer[line * width_band + col] = surface_temperature_vector;
    this->surface_temperature_vector[line][col] = surface_temperature_vector;
  }
};

void Products::short_wave_radiation_function(Reader tal_reader, MTL mtl, uint32_t width_band, int line)
{
  double costheta = sin(mtl.sun_elevation * PI / 180);

  for (int col = 0; col < width_band; col++)
  {
    this->short_wave_radiation_vector[line][col] = (1367 * costheta * tal_reader.read_tiff_pixel(col)) /
                                                   (mtl.distance_earth_sun * mtl.distance_earth_sun);
  }
};

void Products::large_wave_radiation_surface_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    double temperature_pixel = this->surface_temperature_vector[line][col];
    double surface_temperature_pow_4 = temperature_pixel * temperature_pixel * temperature_pixel * temperature_pixel;
    this->large_wave_radiation_surface_vector[line][col] = this->eo_emissivity_vector[line][col] * 5.67 * 1e-8 * surface_temperature_pow_4;
  }
};

void Products::large_wave_radiation_atmosphere_function(uint32_t width_band, double temperature, int line)
{
  double temperature_kelvin = temperature + 273.15;
  double temperature_kelvin_pow_4 = temperature_kelvin * temperature_kelvin * temperature_kelvin * temperature_kelvin;

  for (int col = 0; col < width_band; col++)
    this->large_wave_radiation_atmosphere_vector[line][col] = this->ea_emissivity_vector[line][col] * 5.67 * 1e-8 * temperature_kelvin_pow_4;
};

void Products::net_radiation_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    this->net_radiation_vector[line][col] = this->short_wave_radiation_vector[line][col] - (this->short_wave_radiation_vector[line][col] * this->albedo_vector[line][col]) +
                                            this->large_wave_radiation_atmosphere_vector[line][col] - this->large_wave_radiation_surface_vector[line][col] -
                                            (1 - this->eo_emissivity_vector[line][col]) * this->large_wave_radiation_atmosphere_vector[line][col];

    if (definitelyLessThan(this->net_radiation_vector[line][col], 0))
      this->net_radiation_vector[line][col] = 0;
  }
};

void Products::soil_heat_flux_function(uint32_t width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    if (essentiallyEqual(this->ndvi_vector[line][col], 0) || definitelyGreaterThan(this->ndvi_vector[line][col], 0))
    {
      double ndvi_pixel_pow_4 = this->ndvi_vector[line][col] * this->ndvi_vector[line][col] * this->ndvi_vector[line][col] * this->ndvi_vector[line][col];
      this->soil_heat_vector[line][col] = (this->surface_temperature_vector[line][col] - 273.15) * (0.0038 + 0.0074 * this->albedo_vector[line][col]) *
                                          (1 - 0.98 * ndvi_pixel_pow_4) * this->net_radiation_vector[line][col];
    }
    else
      this->soil_heat_vector[line][col] = 0.5 * this->net_radiation_vector[line][col];

    if (definitelyLessThan(this->soil_heat_vector[line][col], 0))
      this->soil_heat_vector[line][col] = 0;
  }
};

void Products::d0_fuction(int line)
{
  double CD1 = 20.6;
  double HGHT = 4;

  for (int col = 0; col < width_band; col++)
  {

    double pai = this->pai_vector[line][col];
    double cd1_pai_root = sqrt(CD1 * pai);

    double DISP = HGHT * ((1 - (1 / cd1_pai_root)) + (pow(exp(1.0), -cd1_pai_root) / cd1_pai_root));
    if (pai < 0)
    {
      DISP = 0;
    }

    int pos = line * width_band + col;
    this->d0_pointer[pos] = DISP;
  }
};

void Products::kb_function(double ndvi_max, double ndvi_min, int line)
{
  double HGHT = 4;

  double visc = 0.00001461;
  double pr = 0.71;
  double c1 = 0.320;
  double c2 = 0.264;
  double c3 = 15.1;
  double cd = 0.2;
  double ct = 0.01;
  double sf_c = 0.3;
  double sf_d = 2.5;
  double sf_e = 4.0;

  for (int col = 0; col < width_band; col++)
  {
    int pos = line * width_band + col;

    double zom = this->zom_pointer[pos];
    double u_ast_ini_terra = this->ustar_pointer[pos];
    double PAI = this->pai_vector[line][col];

    double Re_star = (u_ast_ini_terra * 0.009) / visc;
    double Ct_star = pow(pr, -0.667) * pow(Re_star, -0.5);
    double beta = c1 - c2 * (exp((cd * -c3 * PAI)));
    double nec_terra = (cd * PAI) / (beta * beta * 2);

    double kb1_fst_part = (cd * VON_KARMAN) / (4 * ct * beta * (1 - exp(nec_terra * -0.5)));
    double kb1_sec_part = (beta * VON_KARMAN * (zom / HGHT)) / Ct_star;
    double kb1s = (pow(Re_star, 0.25) * 2.46) - 2;

    double fc = 1 - pow((this->ndvi_vector[line][col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);
    double fs = 1 - fc;

    double soil_moisture_day_rel = 0.33;

    double SF = sf_c + (1 / (1 + pow(exp(1.0), (sf_d - (sf_e * soil_moisture_day_rel)))));

    this->kb1_pointer[pos] = ((kb1_fst_part * pow(fc, 2)) +
                              (kb1_sec_part * pow(fc, 2) * pow(fs, 2)) +
                              (pow(fs, 2) * kb1s)) *
                             SF;
  }
};

void Products::zom_fuction(double A_ZOM, double B_ZOM, int line)
{
  double HGHT = 4;
  double CD = 0.01;
  double CR = 0.35;
  double PSICORR = 0.2;

  double gama;
  for (int col = 0; col < width_band; col++)
  {
    gama = pow((CD + CR * (this->pai_vector[line][col] / 2)), -0.5);

    if (gama < 3.3)
      gama = 3.3;

    int pos = line * width_band + col;
    this->zom_pointer[pos] = (HGHT - this->d0_pointer[pos]) * pow(exp(1.0), (-VON_KARMAN * gama) + PSICORR);
  }
};

void Products::zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, int width_band, vector<double> &zom_line)
{
  for (int col = 0; col < width_band; col++)
    zom_line[col] = exp(A_ZOM + B_ZOM * ndvi_line[col]);
};

void Products::ustar_fuction(double u10, int line)
{
  double zu = 10;

  for (int col = 0; col < width_band; col++)
  {
    int pos = line * width_band + col;
    double DISP = this->d0_pointer[pos];
    double zom = this->zom_pointer[pos];
    this->ustar_pointer[pos] = (u10 * VON_KARMAN) / std::log((zu - DISP) / zom);
  }
};

void Products::ustar_fuction(double u200, vector<double> zom_line, int width_band, vector<double> &ustar_line)
{
  for (int col = 0; col < width_band; col++)
    ustar_line[col] = (VON_KARMAN * u200) / log(200 / zom_line[col]);
};

void Products::aerodynamic_resistance_fuction(int line)
{
  double zu = 10.0;

  for (int col = 0; col < width_band; col++)
  {
    int pos = line * width_band + col;
    double DISP = this->d0_pointer[pos];
    double zom = this->zom_pointer[pos];
    double zoh_terra = zom / pow(exp(1.0), (this->kb1_pointer[pos]));

    double temp_kb_1_terra = log(zom / zoh_terra);
    double temp_rah1_terra = (1 / (this->ustar_pointer[pos] * VON_KARMAN));
    double temp_rah2 = log(((zu - DISP) / zom));
    double temp_rah3_terra = temp_rah1_terra * temp_kb_1_terra;

    this->aerodynamic_resistance_pointer[pos] = temp_rah1_terra * temp_rah2 + temp_rah3_terra;
  }
};

void Products::aerodynamic_resistance_fuction(vector<double> ustar_line, int width_band, vector<double> &aerodynamic_resistance_line)
{
  for (int col = 0; col < width_band; col++)
    aerodynamic_resistance_line[col] = log(20) / (ustar_line[col] * VON_KARMAN);
};

void Products::sensible_heat_flux_function(double a, double b, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    int pos = line * width_band + col;
    sensible_heat_flux_vector[line][col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_vector[line][col] - 273.15)) / aerodynamic_resistance_pointer[pos];

    if (!isnan(sensible_heat_flux_vector[line][col]) && sensible_heat_flux_vector[line][col] > (net_radiation_vector[line][col] - soil_heat_vector[line][col]))
    {
      sensible_heat_flux_vector[line][col] = net_radiation_vector[line][col] - soil_heat_vector[line][col];
    }
  }
};

void Products::latent_heat_flux_function(int width_band, int line)
{
  for (int col = 0; col < width_band; col++)
    this->latent_heat_flux_vector[line][col] = this->net_radiation_vector[line][col] - this->soil_heat_vector[line][col] - this->sensible_heat_flux_vector[line][col];
};

void Products::net_radiation_24h_function(double Ra24h, double Rs24h, int width_band, int line)
{
  int FL = 110;

  for (int col = 0; col < width_band; col++)
    this->net_radiation_24h_vector[line][col] = (1 - this->albedo_vector[line][col]) * Rs24h - FL * Rs24h / Ra24h;
};

void Products::evapotranspiration_fraction_fuction(int width_band, int line)
{
  for (int col = 0; col < width_band; col++)
  {
    double latent = this->latent_heat_flux_vector[line][col];
    double net_rad = this->net_radiation_vector[line][col];
    double soil_heat = this->soil_heat_vector[line][col];

    this->evapotranspiration_fraction_vector[line][col] = latent / (net_rad - soil_heat);
  }
};

void Products::sensible_heat_flux_24h_fuction(int width_band, int line)
{
  for (int col = 0; col < width_band; col++)
    this->sensible_heat_flux_24h_vector[line][col] = (1 - this->evapotranspiration_fraction_vector[line][col]) * this->net_radiation_24h_vector[line][col];
};

void Products::latent_heat_flux_24h_function(int width_band, int line)
{
  for (int col = 0; col < width_band; col++)
    this->latent_heat_flux_24h_vector[line][col] = this->evapotranspiration_fraction_vector[line][col] * this->net_radiation_24h_vector[line][col];
};

void Products::evapotranspiration_24h_function(Station station, int width_band, int line)
{
  for (int col = 0; col < width_band; col++)
    this->evapotranspiration_24h_vector[line][col] = (this->latent_heat_flux_24h_vector[line][col] * 86400) / ((2.501 - 0.00236 * (station.v7_max + station.v7_min) / 2) * 1e+6);
};

void Products::evapotranspiration_function(int width_band, int line)
{
  for (int col = 0; col < width_band; col++)
    this->evapotranspiration_vector[line][col] = this->net_radiation_24h_vector[line][col] * this->evapotranspiration_fraction_vector[line][col] * 0.035;
};

string Products::rah_correction_function_blocks(double ndvi_min, double ndvi_max, Candidate hot_pixel, Candidate cold_pixel)
{
  system_clock::time_point begin_core, end_core;
  int64_t general_time_core, initial_time_core, final_time_core;

  // ========= CUDA Setup
  int dev = 0;
  cudaDeviceProp deviceProp;
  HANDLE_ERROR(cudaGetDeviceProperties(&deviceProp, dev));
  HANDLE_ERROR(cudaSetDevice(dev));

  int num_threads = 1024;
  int num_blocks = ceil(width_band * height_band / num_threads);

  double hot_pixel_aerodynamic = aerodynamic_resistance_pointer[hot_pixel.line * width_band + hot_pixel.col];
  hot_pixel.aerodynamic_resistance.push_back(hot_pixel_aerodynamic);

  double cold_pixel_aerodynamic = aerodynamic_resistance_pointer[cold_pixel.line * width_band + cold_pixel.col];
  cold_pixel.aerodynamic_resistance.push_back(cold_pixel_aerodynamic);

  double fc_hot = 1 - pow((ndvi_vector[hot_pixel.line][hot_pixel.col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);
  double fc_cold = 1 - pow((ndvi_vector[cold_pixel.line][cold_pixel.col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);

  double *devZom, *devTS, *devUstarR, *devUstarW, *devRahR, *devRahW, *devD0, *devKB1, *devH;
  HANDLE_ERROR(cudaMalloc((void **)&devZom, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devD0, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devKB1, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devTS, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devUstarR, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devUstarW, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devRahR, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devRahW, nBytes_band));
  HANDLE_ERROR(cudaMalloc((void **)&devH, nBytes_band));

  for (int i = 0; i < 2; i++)
  {
    this->rah_ini_pq_terra = hot_pixel.aerodynamic_resistance[i];
    this->rah_ini_pf_terra = cold_pixel.aerodynamic_resistance[i];

    double LEc_terra = 0.55 * fc_hot * (hot_pixel.net_radiation - hot_pixel.soil_heat_flux) * 0.78;
    double LEc_terra_pf = 1.75 * fc_cold * (cold_pixel.net_radiation - cold_pixel.soil_heat_flux) * 0.78;

    this->H_pf_terra = cold_pixel.net_radiation - cold_pixel.soil_heat_flux - LEc_terra_pf;
    double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);

    this->H_pq_terra = hot_pixel.net_radiation - hot_pixel.soil_heat_flux - LEc_terra;
    double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);

    double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
    double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

    HANDLE_ERROR(cudaMemcpy(devTS, surface_temperature_pointer, nBytes_band, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(devZom, zom_pointer, nBytes_band, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(devD0, d0_pointer, nBytes_band, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(devKB1, kb1_pointer, nBytes_band, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(devUstarR, ustar_pointer, nBytes_band, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(devRahR, aerodynamic_resistance_pointer, nBytes_band, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(devH, sensible_heat_flux_pointer, nBytes_band, cudaMemcpyHostToDevice));

    // ==== Paralelization core
    begin_core = system_clock::now();
    initial_time_core = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

    rah_correction_cycle_STEEP<<<num_blocks, num_threads>>>(devTS, devD0, devKB1, devZom, devUstarR, devUstarW, devRahR, devRahW, devH, a, b, height_band, width_band);
    HANDLE_ERROR(cudaDeviceSynchronize());
    HANDLE_ERROR(cudaGetLastError());

    end_core = system_clock::now();
    general_time_core = duration_cast<milliseconds>(end_core - begin_core).count();
    final_time_core = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    // ====

    HANDLE_ERROR(cudaMemcpy(ustar_pointer, devUstarW, nBytes_band, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(aerodynamic_resistance_pointer, devRahW, nBytes_band, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(sensible_heat_flux_pointer, devH, nBytes_band, cudaMemcpyDeviceToHost));

    double rah_hot = this->aerodynamic_resistance_pointer[hot_pixel.line * width_band + hot_pixel.col];
    hot_pixel.aerodynamic_resistance.push_back(rah_hot);

    double rah_cold = this->aerodynamic_resistance_pointer[cold_pixel.line * width_band + cold_pixel.col];
    cold_pixel.aerodynamic_resistance.push_back(rah_cold);
  }

  HANDLE_ERROR(cudaFree(devZom));
  HANDLE_ERROR(cudaFree(devD0));
  HANDLE_ERROR(cudaFree(devKB1));
  HANDLE_ERROR(cudaFree(devTS));
  HANDLE_ERROR(cudaFree(devUstarR));
  HANDLE_ERROR(cudaFree(devUstarW));
  HANDLE_ERROR(cudaFree(devRahR));
  HANDLE_ERROR(cudaFree(devRahW));
  HANDLE_ERROR(cudaDeviceReset());

  return "P2 - RAH - PARALLEL - CORE, " + to_string(general_time_core) + ", " + to_string(initial_time_core) + ", " + to_string(final_time_core) + "\n";
}