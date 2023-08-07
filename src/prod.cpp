#include "prod.h"

/**
 * @brief Calculates the radiance of a line of pixels in a TIFF image.
 *
 * @param bands_resampled An array of pointers to TIFF images, one for each band.
 * @param width_band The width of the TIFF image.
 * @param sample_bands The number of bands in the TIFF image.
 * @param mtl A MTL structure containing metadata about the TIFF image.
 * @param sensor A Sensor structure containing information about the sensor used to capture the TIFF image.
 * @param line The line number to calculate the radiance for.
 * @param radiance_line A 2D array to store the calculated radiance values for each pixel in the line.
 */
void radiance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, double radiance_line[][8])
{
  int inital_index = (mtl.number_sensor == 8) ? 7 : 1;

  for (int i_band = inital_index; i_band < 8; i_band++)
  {
    TIFF *curr_band = bands_resampled[i_band];
    tdata_t line_buff = _TIFFmalloc(TIFFScanlineSize(curr_band));
    PixelReader pixel_reader = read_line_tiff(curr_band, line_buff, line);

    for (int col = 0; col < width_band; col++)
    {
      double band_pixel = pixel_reader.read_pixel(col);

      // TODO: deixar apenas o else (?)
      double radiance_pixel = NaN;
      if (mtl.number_sensor == 8)
        radiance_pixel = max(band_pixel * mtl.rad_mult_10 + mtl.rad_add_10, 0.0);
      else
        radiance_pixel = max(band_pixel * sensor.parameters[i_band][sensor.GRESCALE] + sensor.parameters[i_band][sensor.BRESCALE], 0.0);

      radiance_line[col][i_band] = radiance_pixel;
    }

    _TIFFfree(line_buff);
  }
}

/**
 * @brief Calculates the reflectance of a line of pixels in a TIFF image.
 *
 * @param bands_resampled An array of pointers to TIFF images, one for each band.
 * @param width_band The width of the TIFF image.
 * @param sample_bands The number of bands in the TIFF image.
 * @param mtl A MTL structure containing metadata about the TIFF image.
 * @param sensor A Sensor structure containing information about the sensor used to capture the TIFF image.
 * @param line The line number to calculate the reflectance for.
 * @param reflectance_line A 2D array to store the calculated reflectance values for each pixel in the line.
 */
void reflectance_function(TIFF **bands_resampled, uint32 width_band, uint16 sample_bands, MTL mtl, Sensor sensor, int line, double reflectance_line[][8])
{
  for (int i_band = 1; i_band < 8; i_band++)
  {
    TIFF *curr_band = bands_resampled[i_band];
    tdata_t line_buff = _TIFFmalloc(TIFFScanlineSize(curr_band));
    PixelReader pixel_reader = read_line_tiff(curr_band, line_buff, line);

    for (int col = 0; col < width_band; col++)
    {
      double band_pixel = pixel_reader.read_pixel(col);
      double radiance_pixel = band_pixel * sensor.parameters[i_band][sensor.GRESCALE] + sensor.parameters[i_band][sensor.BRESCALE];

      double reflectance_pixel = NaN;
      if (mtl.number_sensor == 8)
        reflectance_pixel = radiance_pixel / sin(mtl.sun_elevation * PI / 180);
      else
        reflectance_pixel = (PI * radiance_pixel * mtl.distance_earth_sun * mtl.distance_earth_sun) / (sensor.parameters[i_band][sensor.ESUN] * sin(mtl.sun_elevation * PI / 180));

      reflectance_line[col][i_band] = reflectance_pixel;
    }

    _TIFFfree(line_buff);
  }
}

/**
 * @brief Calculates the albedo of a line of pixels in a TIFF image.
 *
 * @param tal_reader A PixelReader object for reading pixel data from the TAL TIFF image.
 * @param sensor A Sensor structure containing information about the sensor used to capture the TIFF image.
 * @param width_band The width of the TIFF image.
 * @param number_sensor The number of bands in the TIFF image.
 * @param reflectance_line A 2D array containing the calculated reflectance values for each pixel in the line.
 * @param albedo_line An array to store the calculated albedo values for each pixel in the line.
 */
void albedo_function(PixelReader tal_reader, Sensor sensor, uint32 width_band, int number_sensor, double reflectance_line[][8], double albedo_line[])
{
  int final_tif_calc = number_sensor == 8 ? 6 : 7;

  for (int col = 0; col < width_band; col++)
  {
    double alb = reflectance_line[col][1] * sensor.parameters[1][sensor.WB] +
                 reflectance_line[col][2] * sensor.parameters[2][sensor.WB] +
                 reflectance_line[col][3] * sensor.parameters[3][sensor.WB] +
                 reflectance_line[col][4] * sensor.parameters[4][sensor.WB] +
                 reflectance_line[col][5] * sensor.parameters[5][sensor.WB] +
                 reflectance_line[col][number_sensor] * sensor.parameters[number_sensor][sensor.WB];

    alb = (alb - 0.03) / (tal_reader.read_pixel(col) * tal_reader.read_pixel(col));

    albedo_line[col] = alb;
  }
}

/**
 * @brief  Normalized Difference Vegetation Index (NDVI) is computed.
 * @note   Values for NDVI range between -1 and 1.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  ndvi_line[]: Auxiliary array for save the calculated value of NDVI for the line.
 */
void ndvi_function(double reflectance_line[][8], uint32 width_band, double ndvi_line[])
{
  for (int col = 0; col < width_band; col++)
  {
    ndvi_line[col] = (reflectance_line[col][4] - reflectance_line[col][3]) /
                     (reflectance_line[col][4] + reflectance_line[col][3]);
  }
};

/**
 * @brief  Plant Area Index (PAI) is computed.
 * @note   Values for PAI range between -1 and 1.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  pai_line[]: Auxiliary array for save the calculated value of NDVI for the line.
 */
void pai_function(double reflectance_line[][8], uint32 width_band, double pai_line[])
{
  for (int col = 0; col < width_band; col++)
  {
    double pai_value = 10.1 * (reflectance_line[col][4] - sqrt(reflectance_line[col][3])) + 3.1;

    if (pai_value < 0)
      pai_value = 0;

    pai_line[col] = pai_value;
  }
};

/**
 * @brief  Leaf Area Index (LAI) is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  lai_line[]: Auxiliary array for save the calculated value of LAI for the line.
 */
void lai_function(double reflectance_line[][8], uint32 width_band, double lai_line[])
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

/**
 * @brief  Enhanced Vegetation Index (EVI) is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  evi_line[]: Auxiliary array for save the calculated value of EVI for the line.
 */
void evi_function(double reflectance_line[][8], uint32 width_band, double evi_line[])
{
  for (int col = 0; col < width_band; col++)
  {
    evi_line[col] = 2.5 * ((reflectance_line[col][4] - reflectance_line[col][3]) /
                           (reflectance_line[col][4] + 6 * reflectance_line[col][3] - 7.5 * reflectance_line[col][1] + 1));
  }
};

/**
 * @brief  Calculates emissivity representing surface behavior for thermal emission in the relatively narrow band 6 of Landsat (10.4 to 12.5 µm),
           expressed as enb.
 * @param  lai_line[]: Array containing the specified line from the LAI computation.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  enb_emissivity_line[]: Auxiliary array for save the calculated value of Enb for the line.
 */
void enb_emissivity_function(double lai_line[], double ndvi_line[], uint32 width_band, double enb_emissivity_line[])
{
  for (int col = 0; col < width_band; col++)
  {
    if (definitelyLessThan(ndvi_line[col], 0) || definitelyGreaterThan(lai_line[col], 2.99))
      enb_emissivity_line[col] = 0.98;
    else
      enb_emissivity_line[col] = 0.97 + 0.0033 * lai_line[col];
  }
};

/**
 * @brief  Calculates emissivity representing surface behavior for thermal emission in the broad thermal spectrum (6 to 14 µm), expressed as eο.
 * @param  lai_line[]: Array containing the specified line from the LAI computation.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  eo_emissivity_line[]: Auxiliary array for save the calculated value of Eo for the line.
 */
void eo_emissivity_function(double lai_line[], double ndvi_line[], uint32 width_band, double eo_emissivity_line[])
{
  for (int col = 0; col < width_band; col++)
  {
    if (definitelyLessThan(ndvi_line[col], 0) || definitelyGreaterThan(lai_line[col], 2.99))
      eo_emissivity_line[col] = 0.98;
    else
      eo_emissivity_line[col] = 0.95 + 0.01 * lai_line[col];
  }
};

/**
 * @brief  Calculates the atmospheric emissivity (ea).
 * @param  tal_line: Object containing the specified line from the tal computation.
 * @param  width_band: Band width.
 * @param  ea_emissivity_line[]: Auxiliary array for save the calculated value of Ea for the line.
 */
void ea_emissivity_function(PixelReader tal_reader, uint32 width_band, double ea_emissivity_line[])
{
  for (int col = 0; col < width_band; col++)
    ea_emissivity_line[col] = 0.85 * pow((-1 * log(tal_reader.read_pixel(col))), 0.09);
};

/**
 * @brief  The surface temperature (TS) is computed.
 * @param  radiance_line[][8]: Radiance for the specific line for each band.
 * @param  enb_emissivity_line[]: Array containing the specified line from the Enb computation.
 * @param  number_sensor: Number of the satellite sensor.
 * @param  width_band: Band width.
 * @param  surface_temperature_line[]: Auxiliary array for save the calculated value of TS for the line.
 */
void surface_temperature_function(double radiance_line[][8], double enb_emissivity_line[], int number_sensor, uint32 width_band, double surface_temperature_line[])
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

/**
 * @brief  Computes Short Wave Radiation (Rs).
 * @param  tal_line: Object containing the specified line from the tal computation.
 * @param  mtl: MTL Struct.
 * @param  width_band: Band width.
 * @param  short_wave_radiation_line[]: Auxiliary array for save the calculated value of Rs for the line.
 */
void short_wave_radiation_function(PixelReader tal_reader, MTL mtl, uint32 width_band, double short_wave_radiation_line[])
{
  double costheta = sin(mtl.sun_elevation * PI / 180);

  for (int col = 0; col < width_band; col++)
  {
    short_wave_radiation_line[col] = (1367 * costheta * tal_reader.read_pixel(col)) /
                                     (mtl.distance_earth_sun * mtl.distance_earth_sun);
  }
};

/**
 * @brief  Computes Large Wave Radiation from Surface (RLSup)
 * @param  eo_emissivity_line[]: Array containing the specified line from the Eo computation.
 * @param  surface_temperature_line[]: Array containing the specified line from the TS computation.
 * @param  width_band: Band width.
 * @param  large_wave_radiation_surface_line[]: Auxiliary array for save the calculated value of RLSup for the line.
 */
void large_wave_radiation_surface_function(double eo_emissivity_line[], double surface_temperature_line[], uint32 width_band, double large_wave_radiation_surface_line[])
{
  for (int col = 0; col < width_band; col++)
  {
    double temperature_pixel = surface_temperature_line[col];
    double surface_temperature_pow_4 = temperature_pixel * temperature_pixel * temperature_pixel * temperature_pixel;
    large_wave_radiation_surface_line[col] = eo_emissivity_line[col] * 5.67 * 1e-8 * surface_temperature_pow_4;
  }
};

/**
 * @brief  Computes Large Wave Radiation from Atmosphere (RLatm)
 * @param  ea_emissivity_line[]: Array containing the specified line from the Ea computation.
 * @param  width_band: Band width.
 * @param  temperature: Near surface air temperature in Kelvin.
 * @param  large_wave_radiation_atmosphere_line[]: Auxiliary array for save the calculated value of RLatm for the line.
 */
void large_wave_radiation_atmosphere_function(double ea_emissivity_line[], uint32 width_band, double temperature, double large_wave_radiation_atmosphere_line[])
{
  double temperature_kelvin = temperature + 273.15;
  double temperature_kelvin_pow_4 = temperature_kelvin * temperature_kelvin * temperature_kelvin * temperature_kelvin;

  for (int col = 0; col < width_band; col++)
    large_wave_radiation_atmosphere_line[col] = ea_emissivity_line[col] * 5.67 * 1e-8 * temperature_kelvin_pow_4;
};

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
                            double eo_emissivity_line[], uint32 width_band, double net_radiation_line[])
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

/**
 * @brief  Computes the Soil heat flux (G).
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  surface_temperature_line[]: Array containing the specified line from the TS computation.
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  width_band: Band width.
 * @param  soil_heat_flux[]: Auxiliary array for save the calculated value of G for the line.
 */
void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], uint32 width_band, double soil_heat_line[])
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
