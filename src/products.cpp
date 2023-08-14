#include "products.h"

void radiance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, vector<vector<double>> &radiance_line)
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
      double band_pixel = pixel_reader.read_tiff_pixel(col);

      // TODO: deixar apenas o else (?)
      double radiance_pixel = NaN;
      if (mtl.number_sensor == 8)
        radiance_pixel = max(band_pixel * mtl.rad_mult_10 + mtl.rad_add_10, 0.0);
      else
        radiance_pixel = max(band_pixel * sensor.parameters[i_band][sensor.GRESCALE] + sensor.parameters[i_band][sensor.BRESCALE], 0.0);

      radiance_line[col][i_band] = radiance_pixel;
    }

    _TIFFfree(band_line_buff);
  }
}

void reflectance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, vector<vector<double>> &reflectance_line)
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
      double band_pixel = pixel_reader.read_tiff_pixel(col);
      double radiance_pixel = band_pixel * sensor.parameters[i_band][sensor.GRESCALE] + sensor.parameters[i_band][sensor.BRESCALE];

      double reflectance_pixel = NaN;
      if (mtl.number_sensor == 8)
        reflectance_pixel = radiance_pixel / sin(mtl.sun_elevation * PI / 180);
      else
        reflectance_pixel = (PI * radiance_pixel * mtl.distance_earth_sun * mtl.distance_earth_sun) / (sensor.parameters[i_band][sensor.ESUN] * sin(mtl.sun_elevation * PI / 180));

      reflectance_line[col][i_band] = reflectance_pixel;
    }

    _TIFFfree(band_line_buff);
  }
}

void albedo_function(Reader tal_reader, Sensor sensor, uint32 width_band, int number_sensor, vector<vector<double>> reflectance_line, vector<double> &albedo_line)
{
  int final_tif_calc = number_sensor == 8 ? 6 : 7;

  for (int col = 0; col < width_band; col++)
  {
    double alb = reflectance_line[col][1] * sensor.parameters[1][sensor.WB] +
                 reflectance_line[col][2] * sensor.parameters[2][sensor.WB] +
                 reflectance_line[col][3] * sensor.parameters[3][sensor.WB] +
                 reflectance_line[col][4] * sensor.parameters[4][sensor.WB] +
                 reflectance_line[col][5] * sensor.parameters[5][sensor.WB] +
                 reflectance_line[col][final_tif_calc] * sensor.parameters[final_tif_calc][sensor.WB];

    alb = (alb - 0.03) / (tal_reader.read_tiff_pixel(col) * tal_reader.read_tiff_pixel(col));

    albedo_line[col] = alb;
  }
}

void ndvi_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &ndvi_line)
{
  for (int col = 0; col < width_band; col++)
  {
    ndvi_line[col] = (reflectance_line[col][4] - reflectance_line[col][3]) /
                     (reflectance_line[col][4] + reflectance_line[col][3]);
  }
};

void pai_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &pai_line)
{
  for (int col = 0; col < width_band; col++)
  {
    double pai_value = 10.1 * (reflectance_line[col][4] - sqrt(reflectance_line[col][3])) + 3.1;

    if (pai_value < 0)
      pai_value = 0;

    pai_line[col] = pai_value;
  }
};

void lai_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &lai_line)
{
  double savi_line[width_band];
  double L = 0.05;

  for (int col = 0; col < width_band; col++)
  {
    savi_line[col] = ((1 + L) * (reflectance_line[col][4] - reflectance_line[col][3])) /
                     (L + (reflectance_line[col][4] + reflectance_line[col][3]));

    if (!isnan(savi_line[col]) && definitelyGreaterThan(savi_line[col], 0.687))
      lai_line[col] = 6;
    else if (!isnan(savi_line[col]) && definitelyLessThan(savi_line[col], 0.1))
      lai_line[col] = 0;
    else
      lai_line[col] = -log((0.69 - savi_line[col]) / 0.59) / 0.91;
  }
};

void evi_function(vector<vector<double>> reflectance_line, uint32 width_band, vector<double> &evi_line)
{
  for (int col = 0; col < width_band; col++)
  {
    evi_line[col] = 2.5 * ((reflectance_line[col][4] - reflectance_line[col][3]) /
                           (reflectance_line[col][4] + 6 * reflectance_line[col][3] - 7.5 * reflectance_line[col][1] + 1));
  }
};

void enb_emissivity_function(vector<double> lai_line, vector<double> ndvi_line, uint32 width_band, vector<double> &enb_emissivity_line)
{
  for (int col = 0; col < width_band; col++)
  {
    if (definitelyLessThan(ndvi_line[col], 0) || definitelyGreaterThan(lai_line[col], 2.99))
      enb_emissivity_line[col] = 0.98;
    else
      enb_emissivity_line[col] = 0.97 + 0.0033 * lai_line[col];
  }
};

void eo_emissivity_function(vector<double> lai_line, vector<double> ndvi_line, uint32 width_band, vector<double> &eo_emissivity_line)
{
  for (int col = 0; col < width_band; col++)
  {
    if (definitelyLessThan(ndvi_line[col], 0) || definitelyGreaterThan(lai_line[col], 2.99))
      eo_emissivity_line[col] = 0.98;
    else
      eo_emissivity_line[col] = 0.95 + 0.01 * lai_line[col];
  }
};

void ea_emissivity_function(Reader tal_reader, uint32 width_band, vector<double> &ea_emissivity_line)
{
  for (int col = 0; col < width_band; col++)
    ea_emissivity_line[col] = 0.85 * pow((-1 * log(tal_reader.read_tiff_pixel(col))), 0.09);
};

void surface_temperature_function(vector<vector<double>> radiance_line, vector<double> enb_emissivity_line, int number_sensor, uint32 width_band, vector<double> &surface_temperature_line)
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
    surface_temperature_line[col] = k2 / (log((enb_emissivity_line[col] * k1 / radiance_line[col][radiance_number]) + 1));
};

void short_wave_radiation_function(Reader tal_reader, MTL mtl, uint32 width_band, vector<double> &short_wave_radiation_line)
{
  double costheta = sin(mtl.sun_elevation * PI / 180);

  for (int col = 0; col < width_band; col++)
  {
    short_wave_radiation_line[col] = (1367 * costheta * tal_reader.read_tiff_pixel(col)) /
                                     (mtl.distance_earth_sun * mtl.distance_earth_sun);
  }
};

void large_wave_radiation_surface_function(vector<double> eo_emissivity_line, vector<double> surface_temperature_line, uint32 width_band, vector<double> &large_wave_radiation_surface_line)
{
  for (int col = 0; col < width_band; col++)
  {
    double temperature_pixel = surface_temperature_line[col];
    double surface_temperature_pow_4 = temperature_pixel * temperature_pixel * temperature_pixel * temperature_pixel;
    large_wave_radiation_surface_line[col] = eo_emissivity_line[col] * 5.67 * 1e-8 * surface_temperature_pow_4;
  }
};

void large_wave_radiation_atmosphere_function(vector<double> ea_emissivity_line, uint32 width_band, double temperature, vector<double> &large_wave_radiation_atmosphere_line)
{
  double temperature_kelvin = temperature + 273.15;
  double temperature_kelvin_pow_4 = temperature_kelvin * temperature_kelvin * temperature_kelvin * temperature_kelvin;

  for (int col = 0; col < width_band; col++)
    large_wave_radiation_atmosphere_line[col] = ea_emissivity_line[col] * 5.67 * 1e-8 * temperature_kelvin_pow_4;
};

void net_radiation_function(vector<double> short_wave_radiation_line, vector<double> large_wave_radiation_surface_line,
                            vector<double> large_wave_radiation_atmosphere_line, vector<double> albedo_line,
                            vector<double> eo_emissivity_line, uint32 width_band, vector<double> &net_radiation_line)
{
  for (int col = 0; col < width_band; col++)
  {
    net_radiation_line[col] = short_wave_radiation_line[col] - (short_wave_radiation_line[col] * albedo_line[col]) +
                              large_wave_radiation_atmosphere_line[col] - large_wave_radiation_surface_line[col] -
                              (1 - eo_emissivity_line[col]) * large_wave_radiation_atmosphere_line[col];

    if (definitelyLessThan(net_radiation_line[col], 0))
      net_radiation_line[col] = 0;
  }
};

void soil_heat_flux_function(vector<double> ndvi_line, vector<double> surface_temperature_line, vector<double> albedo_line, vector<double> net_radiation_line, uint32 width_band, vector<double> &soil_heat_line)
{
  for (int col = 0; col < width_band; col++)
  {
    if (essentiallyEqual(ndvi_line[col], 0) || definitelyGreaterThan(ndvi_line[col], 0))
    {
      double ndvi_pixel_pow_4 = ndvi_line[col] * ndvi_line[col] * ndvi_line[col] * ndvi_line[col];
      soil_heat_line[col] = (surface_temperature_line[col] - 273.15) * (0.0038 + 0.0074 * albedo_line[col]) *
                            (1 - 0.98 * ndvi_pixel_pow_4) * net_radiation_line[col];
    }
    else
      soil_heat_line[col] = 0.5 * net_radiation_line[col];

    if (definitelyLessThan(soil_heat_line[col], 0))
      soil_heat_line[col] = 0;
  }
};

void d0_fuction(vector<double> pai_line, int width_band, vector<double> &d0_line)
{
  double CD1 = 20.6;
  double HGHT = 4;

  for (int col = 0; col < width_band; col++)
  {

    double pai = pai_line[col];
    double cd1_pai_root = sqrt(CD1 * pai);

    double DISP = HGHT * ((1 - (1 / cd1_pai_root)) + (pow(exp(1.0), -cd1_pai_root) / cd1_pai_root));
    if (pai < 0)
    {
      DISP = 0;
    }

    d0_line[col] = DISP;
  }
};

void kb_function(vector<double> ustar_line, vector<double> zom_line, vector<double> pai_line, vector<double> ndvi_line, double ndvi_max, double ndvi_min, int width_band, vector<double> &kb1_line)
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
    double zom = zom_line[col];
    double u_ast_ini_terra = ustar_line[col];
    double PAI = pai_line[col];

    double Re_star = (u_ast_ini_terra * 0.009) / visc;
    double Ct_star = pow(pr, -0.667) * pow(Re_star, -0.5);
    double beta = c1 - c2 * (exp((cd * -c3 * PAI)));
    double nec_terra = (cd * PAI) / (beta * beta * 2);

    double kb1_fst_part = (cd * VON_KARMAN) / (4 * ct * beta * (1 - exp(nec_terra * -0.5)));
    double kb1_sec_part = (beta * VON_KARMAN * (zom / HGHT)) / Ct_star;
    double kb1s = (pow(Re_star, 0.25) * 2.46) - 2;

    double fc = 1 - pow((ndvi_line[col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);
    double fs = 1 - fc;

    double soil_moisture_day_rel = 0.33;

    double SF = sf_c + (1 / (1 + pow(exp(1.0), (sf_d - (sf_e * soil_moisture_day_rel)))));

    kb1_line[col] = ((kb1_fst_part * pow(fc, 2)) +
                     (kb1_sec_part * pow(fc, 2) * pow(fs, 2)) +
                     (pow(fs, 2) * kb1s)) *
                    SF;
  }
};

void zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, vector<double> d0_line, vector<double> pai_line, int width_band, vector<double> &zom_line)
{
  double HGHT = 4;
  double CD = 0.01;
  double CR = 0.35;
  double PSICORR = 0.2;

  double gama;
  for (int col = 0; col < width_band; col++)
  {
    gama = pow((CD + CR * (pai_line[col] / 2)), -0.5);

    if (gama < 3.3)
      gama = 3.3;

    zom_line[col] = (HGHT - d0_line[col]) * pow(exp(1.0), (-VON_KARMAN * gama) + PSICORR);
  }
};

void zom_fuction(double A_ZOM, double B_ZOM, vector<double> ndvi_line, int width_band, vector<double> &zom_line)
{
  for (int col = 0; col < width_band; col++)
    zom_line[col] = exp(A_ZOM + B_ZOM * ndvi_line[col]);
};

void ustar_fuction(double u10, vector<double> zom_line, vector<double> d0_line, int width_band, vector<double> &ustar_line)
{
  double zu = 10;

  for (int col = 0; col < width_band; col++)
  {
    double DISP = d0_line[col];
    double zom = zom_line[col];
    ustar_line[col] = (u10 * VON_KARMAN) / std::log((zu - DISP) / zom);
  }
};

void ustar_fuction(double u200, vector<double> zom_line, int width_band, vector<double> &ustar_line)
{
  for (int col = 0; col < width_band; col++)
    ustar_line[col] = (VON_KARMAN * u200) / log(200 / zom_line[col]);
};

void aerodynamic_resistance_fuction(vector<double> ustar_line, vector<double> zom_line, vector<double> d0_line, vector<double> kb1_line, int width_band, vector<double> &aerodynamic_resistance_line)
{
  double zu = 10.0;

  for (int col = 0; col < width_band; col++)
  {
    double DISP = d0_line[col];
    double zom = zom_line[col];
    double zoh_terra = zom / pow(exp(1.0), (kb1_line[col]));

    double temp_kb_1_terra = log(zom / zoh_terra);
    double temp_rah1_terra = (1 / (ustar_line[col] * VON_KARMAN));
    double temp_rah2 = log(((zu - DISP) / zom));
    double temp_rah3_terra = temp_rah1_terra * temp_kb_1_terra;

    aerodynamic_resistance_line[col] = temp_rah1_terra * temp_rah2 + temp_rah3_terra;
  }
};

void aerodynamic_resistance_fuction(vector<double> ustar_line, int width_band, vector<double> &aerodynamic_resistance_line)
{
  for (int col = 0; col < width_band; col++)
    aerodynamic_resistance_line[col] = log(20) / (ustar_line[col] * VON_KARMAN);
};

void sensible_heat_function_STEEP(Candidate hot_pixel, Candidate cold_pixel, Station station, uint32 height_band, uint32 width_band, vector<vector<double>> ndvi_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> pai_vector, vector<vector<double>> &sensible_heat_flux_vector)
{
  // ============== COMPUTE INITIAL RAH

  double ustar_station = (VON_KARMAN * station.v6) / (log(station.WIND_SPEED / station.SURFACE_ROUGHNESS));
  double u10 = (ustar_station / VON_KARMAN) * log(10 / station.SURFACE_ROUGHNESS);

  double ndvi_min = 1.0;
  double ndvi_max = 0.0;
  for (int line = 0; line < height_band; line++)
  {
    vector<double> ndvi_line = ndvi_vector[line];
    for (int col = 0; col < width_band; col++)
    {
      if (ndvi_line[col] < ndvi_min)
        ndvi_min = ndvi_line[col];
      if (ndvi_line[col] > ndvi_max)
        ndvi_max = ndvi_line[col];
    }
  }

  vector<vector<double>> d0_vector(height_band, vector<double>(width_band));
  vector<vector<double>> zom_vector(height_band, vector<double>(width_band));
  vector<vector<double>> ustar_vector(height_band, vector<double>(width_band));
  vector<vector<double>> kb1_vector(height_band, vector<double>(width_band));
  vector<vector<double>> aerodynamic_resistance_vector(height_band, vector<double>(width_band));

  for (int line = 0; line < height_band; line++)
  {
    d0_fuction(pai_vector[line], width_band, d0_vector[line]);
    zom_fuction(station.A_ZOM, station.B_ZOM, ndvi_vector[line], d0_vector[line], pai_vector[line], width_band, zom_vector[line]);
    ustar_fuction(u10, zom_vector[line], d0_vector[line], width_band, ustar_vector[line]);
    kb_function(ustar_vector[line], zom_vector[line], pai_vector[line], ndvi_vector[line], ndvi_max, ndvi_min, width_band, kb1_vector[line]);
    aerodynamic_resistance_fuction(ustar_vector[line], zom_vector[line], d0_vector[line], kb1_vector[line], width_band, aerodynamic_resistance_vector[line]);
  }

  // ============== COMPUTE FINAL RAH

  vector<vector<double>> ustar_previous(height_band, vector<double>(width_band));
  vector<vector<double>> aerodynamic_resistance_previous(height_band, vector<double>(width_band));

  // Auxiliaries arrays calculation
  double L[width_band];
  double y_01_line[width_band], y_2_line[width_band], x_200_line[width_band];
  double psi_01_line[width_band], psi_2_line[width_band], psi_200_line[width_band];

  double hot_pixel_aerodynamic = aerodynamic_resistance_vector[hot_pixel.line][hot_pixel.col];
  hot_pixel.aerodynamic_resistance.push_back(hot_pixel_aerodynamic);

  double cold_pixel_aerodynamic = aerodynamic_resistance_vector[cold_pixel.line][cold_pixel.col];
  cold_pixel.aerodynamic_resistance.push_back(cold_pixel_aerodynamic);

  double fc_hot = 1 - pow((ndvi_vector[hot_pixel.line][hot_pixel.col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);
  double fc_cold = 1 - pow((ndvi_vector[cold_pixel.line][cold_pixel.col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);

  double H_pf_terra;
  double H_pq_terra;
  double rah_ini_pq_terra;
  double rah_ini_pf_terra;

  for (int i = 0; i < 2; i++)
  {
    ustar_previous = ustar_vector;
    aerodynamic_resistance_previous = aerodynamic_resistance_vector;

    rah_ini_pq_terra = hot_pixel.aerodynamic_resistance[i];
    rah_ini_pf_terra = cold_pixel.aerodynamic_resistance[i];

    double LEc_terra = 0.55 * fc_hot * (hot_pixel.net_radiation - hot_pixel.soil_heat_flux) * 0.78;
    double LEc_terra_pf = 1.75 * fc_cold * (cold_pixel.net_radiation - cold_pixel.soil_heat_flux) * 0.78;

    H_pf_terra = cold_pixel.net_radiation - cold_pixel.soil_heat_flux - LEc_terra_pf;
    double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);

    H_pq_terra = hot_pixel.net_radiation - hot_pixel.soil_heat_flux - LEc_terra;
    double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);

    double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
    double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

    for (int line = 0; line < height_band; line++)
    {
      for (int col = 0; col < width_band; col++)
      {
        double DISP = d0_vector[line][col];
        double dT_ini_terra = (a + b * (surface_temperature_vector[line][col] - 273.15));

        // H_ini_terra
        sensible_heat_flux_vector[line][col] = RHO * SPECIFIC_HEAT_AIR * (dT_ini_terra) / aerodynamic_resistance_previous[line][col];

        // L_MB_terra
        double ustar_pow_3 = ustar_previous[line][col] * ustar_previous[line][col] * ustar_previous[line][col];
        L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_vector[line][col]) / (VON_KARMAN * GRAVITY * sensible_heat_flux_vector[line][col]));

        y_01_line[col] = pow((1 - (16 * 0.1) / L[col]), 0.25);
        y_2_line[col] = pow((1 - (16 * (10 - DISP)) / L[col]), 0.25);
        x_200_line[col] = pow((1 - (16 * (10 - DISP)) / L[col]), 0.25);

        // psi_m200_terra
        if (L[col] < 0)
          psi_200_line[col] = 2 * log((1 + x_200_line[col]) / 2) + log((1 + x_200_line[col] * x_200_line[col]) / 2) - 2 * atan(x_200_line[col]) + 0.5 * PI;
        else
          psi_200_line[col] = -5 * ((10 - DISP) / L[col]);

        // psi_m2_terra
        if (L[col] < 0)
          psi_2_line[col] = 2 * log((1 + y_2_line[col] * y_2_line[col]) / 2);
        else
          psi_2_line[col] = -5 * ((10 - DISP) / L[col]);

        // u*
        double ust = (VON_KARMAN * ustar_previous[line][col]) / (log((10 - DISP) / zom_vector[line][col]) - psi_200_line[col]);

        // rah
        double zoh_terra = zom_vector[line][col] / pow(exp(1.0), (kb1_vector[line][col]));
        double temp_rah1_corr_terra = (ust * VON_KARMAN);
        double temp_rah2_corr_terra = log((10 - DISP) / zom_vector[line][col]) - psi_2_line[col];
        double temp_rah3_corr_terra = temp_rah1_corr_terra * log(zom_vector[line][col] / zoh_terra);
        double rah = (temp_rah1_corr_terra * temp_rah2_corr_terra) + temp_rah3_corr_terra;

        ustar_vector[line][col] = ust;
        aerodynamic_resistance_vector[line][col] = rah;

        if (line == hot_pixel.line && col == hot_pixel.col)
        {
          hot_pixel.aerodynamic_resistance.push_back(rah);
        }

        if (line == cold_pixel.line && col == cold_pixel.col)
        {
          cold_pixel.aerodynamic_resistance.push_back(rah);
        }
      }
    }
  }

  // ============== COMPUTE H

  double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);
  double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);

  double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
  double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

  for (int line = 0; line < height_band; line++)
  {
    for (int col = 0; col < width_band; col++)
    {
      sensible_heat_flux_vector[line][col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_vector[line][col] - 273.15)) / aerodynamic_resistance_vector[line][col];

      if (!isnan(sensible_heat_flux_vector[line][col]) && sensible_heat_flux_vector[line][col] > (net_radiation_vector[line][col] - soil_heat_vector[line][col]))
      {
        sensible_heat_flux_vector[line][col] = net_radiation_vector[line][col] - soil_heat_vector[line][col];
      }
    }
  }
};

void sensible_heat_function_default(Candidate hot_pixel, Candidate cold_pixel, Station station, uint32 height_band, uint32 width_band, vector<vector<double>> ndvi_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> &sensible_heat_flux_vector)
{
  // ============== COMPUTE INITIAL RAH

  double ustar_station = (VON_KARMAN * station.v6) / (log(station.WIND_SPEED / station.SURFACE_ROUGHNESS));
  double u200 = (ustar_station / VON_KARMAN) * log(200 / station.SURFACE_ROUGHNESS);

  vector<vector<double>> d0_vector(height_band, vector<double>(width_band));
  vector<vector<double>> zom_vector(height_band, vector<double>(width_band));
  vector<vector<double>> ustar_vector(height_band, vector<double>(width_band));
  vector<vector<double>> kb1_vector(height_band, vector<double>(width_band));
  vector<vector<double>> aerodynamic_resistance_vector(height_band, vector<double>(width_band));

  for (int line = 0; line < height_band; line++)
  {
    zom_fuction(station.A_ZOM, station.B_ZOM, ndvi_vector[line], width_band, zom_vector[line]);
    ustar_fuction(u200, zom_vector[line], width_band, ustar_vector[line]);
    aerodynamic_resistance_fuction(ustar_vector[line], width_band, aerodynamic_resistance_vector[line]);
  }

  // ============== COMPUTE FINAL RAH

  vector<vector<double>> ustar_previous(height_band, vector<double>(width_band));
  vector<vector<double>> aerodynamic_resistance_previous(height_band, vector<double>(width_band));

  // Auxiliaries arrays calculation
  double L[width_band];
  double y_01_line[width_band], y_2_line[width_band], x_200_line[width_band];
  double psi_01_line[width_band], psi_2_line[width_band], psi_200_line[width_band];

  double hot_pixel_aerodynamic = aerodynamic_resistance_vector[hot_pixel.line][hot_pixel.col];
  hot_pixel.aerodynamic_resistance.push_back(hot_pixel_aerodynamic);

  double cold_pixel_aerodynamic = aerodynamic_resistance_vector[cold_pixel.line][cold_pixel.col];
  cold_pixel.aerodynamic_resistance.push_back(cold_pixel_aerodynamic);

  double H_pq_terra;
  double H_pf_terra;
  double rah_ini_pq_terra;
  double rah_ini_pf_terra;
  double rah_final_pq_terra;
  double rah_final_pf_terra;

  int i = 0;
  bool Error = true;
  while (Error)
  {
    ustar_previous = ustar_vector;
    aerodynamic_resistance_previous = aerodynamic_resistance_vector;

    H_pq_terra = hot_pixel.net_radiation - hot_pixel.soil_heat_flux;
    H_pf_terra = cold_pixel.net_radiation - cold_pixel.soil_heat_flux;

    rah_ini_pq_terra = hot_pixel.aerodynamic_resistance[i];
    rah_ini_pf_terra = cold_pixel.aerodynamic_resistance[i];

    double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);
    double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);

    double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
    double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

    for (int line = 0; line < height_band; line++)
    {
      for (int col = 0; col < width_band; col++)
      {
        sensible_heat_flux_vector[line][col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_vector[line][col] - 273.15)) / aerodynamic_resistance_previous[line][col];

        double ustar_pow_3 = ustar_previous[line][col] * ustar_previous[line][col] * ustar_previous[line][col];
        L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_vector[line][col]) / (VON_KARMAN * GRAVITY * sensible_heat_flux_vector[line][col]));

        y_01_line[col] = pow((1 - (16 * 0.1) / L[col]), 0.25);
        y_2_line[col] = pow((1 - (16 * 2) / L[col]), 0.25);
        x_200_line[col] = pow((1 - (16 * 200) / L[col]), 0.25);

        if (!isnan(L[col]) && L[col] > 0)
          psi_01_line[col] = -5 * (0.1 / L[col]);
        else
          psi_01_line[col] = 2 * log((1 + y_01_line[col] * y_01_line[col]) / 2);

        if (!isnan(L[col]) && L[col] > 0)
          psi_2_line[col] = -5 * (2 / L[col]);
        else
          psi_2_line[col] = 2 * log((1 + y_2_line[col] * y_2_line[col]) / 2);

        if (!isnan(L[col]) && L[col] > 0)
          psi_200_line[col] = -5 * (2 / L[col]);
        else
          psi_200_line[col] = 2 * log((1 + x_200_line[col]) / 2) + log((1 + x_200_line[col] * x_200_line[col]) / 2) - 2 * atan(x_200_line[col]) + 0.5 * PI;

        double ust = (VON_KARMAN * u200) / (log(200 / zom_vector[line][col]) - psi_200_line[col]);
        ustar_vector[line][col] = ust;

        rah_final_pq_terra = (log(2 / 0.1) - psi_2_line[col] + psi_01_line[col]) / (ust * VON_KARMAN);
        aerodynamic_resistance_vector[line][col] = rah_final_pq_terra;

        if (line == hot_pixel.line && col == hot_pixel.col)
        {
          hot_pixel.aerodynamic_resistance.push_back(rah_final_pq_terra);
        }
      }
    }

    Error = (fabs(1 - rah_ini_pq_terra / rah_final_pq_terra) >= 0.05);
    i++;
  }

  // ============== COMPUTE H

  double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);
  double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);

  double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
  double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

  for (int line = 0; line < height_band; line++)
  {
    for (int col = 0; col < width_band; col++)
    {
      sensible_heat_flux_vector[line][col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_vector[line][col] - 273.15)) / aerodynamic_resistance_vector[line][col];

      if (!isnan(sensible_heat_flux_vector[line][col]) && sensible_heat_flux_vector[line][col] > (net_radiation_vector[line][col] - soil_heat_vector[line][col]))
      {
        sensible_heat_flux_vector[line][col] = net_radiation_vector[line][col] - soil_heat_vector[line][col];
      }
    }
  }
};

void latent_heat_flux_function(vector<double> net_radiation_line, vector<double> soil_heat_flux_line, vector<double> sensible_heat_flux_line, int width_band, vector<double> &latent_heat_flux)
{
  for (int col = 0; col < width_band; col++)
    latent_heat_flux[col] = net_radiation_line[col] - soil_heat_flux_line[col] - sensible_heat_flux_line[col];
};

void net_radiation_24h_function(vector<double> albedo_line, double Ra24h, double Rs24h, int width_band, vector<double> &net_radiation_24h_line)
{
  int FL = 110;

  for (int col = 0; col < width_band; col++)
    net_radiation_24h_line[col] = (1 - albedo_line[col]) * Rs24h - FL * Rs24h / Ra24h;
};

void evapotranspiration_fraction_fuction(vector<double> latent_heat_flux_line, vector<double> net_radiation_line, vector<double> soil_heat_line, int width_band, vector<double> &evapotranspiration_fraction_line)
{
  for (int col = 0; col < width_band; col++)
  {
    double latent = latent_heat_flux_line[col];
    double net_rad = net_radiation_line[col];
    double soil_heat = soil_heat_line[col];

    evapotranspiration_fraction_line[col] = latent / (net_rad - soil_heat);
  }
};

void sensible_heat_flux_24h_fuction(vector<double> evapotranspiration_fraction_line, vector<double> net_radiation_24h_line, int width_band, vector<double> &sensible_heat_flux_24h_line)
{
  for (int col = 0; col < width_band; col++)
    sensible_heat_flux_24h_line[col] = (1 - evapotranspiration_fraction_line[col]) * net_radiation_24h_line[col];
};

void latent_heat_flux_24h_function(vector<double> evapotranspiration_fraction_line, vector<double> net_radiation_24h_line, int width_band, vector<double> &latent_heat_flux_24h_line)
{
  for (int col = 0; col < width_band; col++)
    latent_heat_flux_24h_line[col] = evapotranspiration_fraction_line[col] * net_radiation_24h_line[col];
};

void evapotranspiration_24h_function(vector<double> latent_heat_flux_24h_line, Station station, int width_band, vector<double> &evapotranspiration_24h_line)
{
  for (int col = 0; col < width_band; col++)
    evapotranspiration_24h_line[col] = (latent_heat_flux_24h_line[col] * 86400) / ((2.501 - 0.00236 * (station.v7_max + station.v7_min) / 2) * 1e+6);
};

void evapotranspiration_function(vector<double> net_radiation_24h_line, vector<double> evapotranspiration_fraction_line, int width_band, vector<double> &evapotranspiration_24h_line)
{
  for (int col = 0; col < width_band; col++)
    evapotranspiration_24h_line[col] = net_radiation_24h_line[col] * evapotranspiration_fraction_line[col] * 0.035;
};
