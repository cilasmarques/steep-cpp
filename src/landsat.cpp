#include "landsat.h"

/**
 * @brief  Empty constructor.
 */
Landsat::Landsat(){};

/**
 * @brief  Constructor of the struct.
 * @param  tal_path: Path to tal TIFF.
 * @param  output_path: Output path where TIFF should be saved.
 * @param  method: Pixel selection method.
 * @param  noData: TIFFs no data value.
 * @param  land_cover_path: Path to Land Cover TIFF.
 */
Landsat::Landsat(string tal_path, string output_path, int method, double noData, string land_cover_path)
{
  this->tal_path = tal_path;
  this->output_path = output_path;
  this->land_cover_path = land_cover_path;

  // Initialize the path of products TIFF based on the output path.
  this->albedo_path = output_path + "/alb.tif";
  this->ndvi_path = output_path + "/NDVI.tif";
  this->evi_path = output_path + "/EVI.tif";
  this->lai_path = output_path + "/LAI.tif";
  this->pai_path = output_path + "/PAI.tif";
  this->soil_heat_path = output_path + "/G.tif";
  this->surface_temperature_path = output_path + "/TS.tif";
  this->net_radiation_path = output_path + "/Rn.tif";
  this->evapotranspiration_fraction_path = output_path + "/EF.tif";
  this->evapotranspiration_24h_path = output_path + "/ET24h.tif";
  this->zom_path = output_path + "/zom.tif";
  this->ustar_path = output_path + "/ustar.tif";
  this->aerodynamic_resistance_path = output_path + "/Rah.tif";
  this->sensible_heat_flux_path = output_path + "/H.tif";
  this->ustar_tif0_path = output_path + "/ustar_tif0.tif";
  this->ustar_tif1_path = output_path + "/ustar_tif1.tif";
  this->aerodynamic_resistance_tif0_path = output_path + "/Rah_tif0.tif";
  this->aerodynamic_resistance_tif1_path = output_path + "/Rah_tif1.tif";
  this->latent_heat_flux_path = output_path + "/LatentHF.tif";
  this->net_radiation_24h_path = output_path + "/Rn24h.tif";
  this->latent_heat_flux_24h_path = output_path + "/LatentHF24h.tif";
  this->noData = noData;
  this->method = method;
};

/**
 * @brief  Initializes TIFFs of the partial execution products as writable. Doing their setup based upon the tal TIFF characteristics.
 * @param  **tal: Tal TIFF.
 * @param  **albedo: Albedo TIFF.
 * @param  **ndvi: NDVI TIFF.
 * @param  **evi: EVI TIFF.
 * @param  **lai: LAI TIFF.
 * @param  **soil_heat: Soil heat flux TIFF.
 * @param  **surface_temperature: Surface temperature TIFF.
 * @param  **net_radiation: Net radiation TIFF.
 */
void Landsat::create_tiffs(TIFF **tal, TIFF **albedo, TIFF **ndvi, TIFF **evi, TIFF **lai, TIFF **pai, TIFF **soil_heat, TIFF **surface_temperature, TIFF **net_radiation)
{
  *albedo = TIFFOpen(albedo_path.c_str(), "w8m");
  setup(*albedo, *tal);

  *ndvi = TIFFOpen(ndvi_path.c_str(), "w8m");
  setup(*ndvi, *tal);

  *evi = TIFFOpen(evi_path.c_str(), "w8m");
  setup(*evi, *tal);

  *lai = TIFFOpen(lai_path.c_str(), "w8m");
  setup(*lai, *tal);

  *pai = TIFFOpen(pai_path.c_str(), "w8m");
  setup(*pai, *tal);

  *soil_heat = TIFFOpen(soil_heat_path.c_str(), "w8m");
  setup(*soil_heat, *tal);

  *surface_temperature = TIFFOpen(surface_temperature_path.c_str(), "w8m");
  setup(*surface_temperature, *tal);

  *net_radiation = TIFFOpen(net_radiation_path.c_str(), "w8m");
  setup(*net_radiation, *tal);
};

/**
 * @brief  Calculates the partials products (e. g. Albedo, NDVI, Rn, G) of the SEBAL execution.
 * @param  read_bands[]: Satellite images as TIFFs.
 * @param  mtl: MTL struct.
 * @param  station: Station struct.
 * @param  sensor: Sensor struct.
 */
void Landsat::process_partial_products(TIFF *bands_resampled[], MTL mtl, Station station, Sensor sensor)
{
  uint16 sample_bands;
  uint32 height_band, width_band;
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGELENGTH, &height_band);
  TIFFGetField(bands_resampled[1], TIFFTAG_IMAGEWIDTH, &width_band);
  TIFFGetField(bands_resampled[1], TIFFTAG_SAMPLEFORMAT, &sample_bands);

  TIFF *tal = TIFFOpen(this->tal_path.c_str(), "rm");
  check_open_tiff(tal);

  TIFF *albedo, *ndvi, *evi, *pai, *lai, *soil_heat, *surface_temperature, *net_radiation;
  create_tiffs(&tal, &albedo, &ndvi, &evi, &pai, &lai, &soil_heat, &surface_temperature, &net_radiation);

  // Declare auxiliaries arrays
  double radiance_line[width_band][8];
  double reflectance_line[width_band][8];

  double albedo_line[width_band];
  double ndvi_line[width_band];
  double lai_line[width_band];
  double evi_line[width_band];
  double pai_line[width_band];
  double enb_emissivity_line[width_band];
  double eo_emissivity_line[width_band];
  double ea_emissivity_line[width_band];
  double short_wave_radiation_line[width_band];
  double surface_temperature_line[width_band];
  double large_wave_radiation_surface_line[width_band];
  double large_wave_radiation_atmosphere_line[width_band];
  double net_radiation_line[width_band];
  double soil_heat_line[width_band];

  for (int line = 0; line < height_band; line++)
  {
    tdata_t tal_line_buff = _TIFFmalloc(TIFFScanlineSize(tal));
    PixelReader tal_reader = read_line_tiff(tal, tal_line_buff, line);

    radiance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line, radiance_line);
    reflectance_function(bands_resampled, width_band, sample_bands, mtl, sensor, line, reflectance_line);
    albedo_function(tal_reader, sensor, width_band, mtl.number_sensor, reflectance_line, albedo_line);

    // Vegetation indices
    ndvi_function(reflectance_line, width_band, ndvi_line);
    pai_function(reflectance_line, width_band, pai_line);
    lai_function(reflectance_line, width_band, lai_line);
    evi_function(reflectance_line, width_band, evi_line);

    // Emissivity indices
    ea_emissivity_function(tal_reader, width_band, ea_emissivity_line);
    enb_emissivity_function(lai_line, ndvi_line, width_band, enb_emissivity_line);
    eo_emissivity_function(lai_line, ndvi_line, width_band, eo_emissivity_line);
    surface_temperature_function(radiance_line, enb_emissivity_line, mtl.number_sensor, width_band, surface_temperature_line);

    // Radiation waves
    short_wave_radiation_function(tal_reader, mtl, width_band, short_wave_radiation_line);
    large_wave_radiation_surface_function(eo_emissivity_line, surface_temperature_line, width_band, large_wave_radiation_surface_line);
    large_wave_radiation_atmosphere_function(ea_emissivity_line, width_band, station.temperature_image, large_wave_radiation_atmosphere_line);

    // Main products
    net_radiation_function(short_wave_radiation_line, large_wave_radiation_surface_line, large_wave_radiation_atmosphere_line, albedo_line, eo_emissivity_line, width_band, net_radiation_line);
    soil_heat_flux_function(ndvi_line, surface_temperature_line, albedo_line, net_radiation_line, width_band, soil_heat_line);

    _TIFFfree(tal_line_buff);

    TIFFWriteScanline(albedo, albedo_line, line);
    TIFFWriteScanline(ndvi, ndvi_line, line);
    TIFFWriteScanline(evi, evi_line, line);
    TIFFWriteScanline(lai, lai_line, line);
    TIFFWriteScanline(pai, pai_line, line);
    TIFFWriteScanline(soil_heat, soil_heat_line, line);
    TIFFWriteScanline(surface_temperature, surface_temperature_line, line);
    TIFFWriteScanline(net_radiation, net_radiation_line, line);
  }

  TIFFClose(albedo);
  TIFFClose(ndvi);
  TIFFClose(evi);
  TIFFClose(lai);
  TIFFClose(pai);
  TIFFClose(soil_heat);
  TIFFClose(surface_temperature);
  TIFFClose(net_radiation);
  TIFFClose(tal);
};

