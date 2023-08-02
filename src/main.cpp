#include <fstream>
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"

int main(int argc, char *argv[])
{
  const int INPUT_BAND_BLUE_INDEX = 1;
  const int INPUT_BAND_GREEN_INDEX = 2;
  const int INPUT_BAND_RED_INDEX = 3;
  const int INPUT_BAND_NIR_INDEX = 4;
  const int INPUT_BAND_SWIR1_INDEX = 5;
  const int INPUT_BAND_TERMAL_INDEX = 6;
  const int INPUT_BAND_SWIR2_INDEX = 7;
  const int INPUT_BAND_TAL_INDEX = 8;
  const int INPUT_MTL_DATA_INDEX = 9;
  const int INPUT_STATION_DATA_INDEX = 10;
  // const int INPUT_LAND_COVER_INDEX = 11;
  // const int OUTPUT_FOLDER = 12;
  // const int INPUT_ALG_METHOD_INDEX = 13;
  // const int INPUT_NAN_VALUE_INDEX = 14;


  // ================= LOAD DATA

  //load meta data
  string path_meta_file = argv[INPUT_MTL_DATA_INDEX];
  MTL mtl = MTL(path_meta_file);

  //load sensor data
  string station_data_path = argv[INPUT_STATION_DATA_INDEX];
  Station station = Station(station_data_path, mtl.image_hour);

  //load sensor data
  Sensor sensor = Sensor(mtl.number_sensor, mtl.year);


  // ================= DEFINE GLOBAL VARIABLES

  double cos_theta = sin(mtl.sun_elevation * PI / 180);

  string path_tiff_band_red = argv[INPUT_BAND_RED_INDEX];
  TIFF *band_red = TIFFOpen(path_tiff_band_red.c_str(), "rm");

  uint32 height_band, width_band;
  uint16 sample_bands;

  TIFFGetField(band_red, TIFFTAG_IMAGELENGTH, &height_band);
  TIFFGetField(band_red, TIFFTAG_IMAGEWIDTH, &width_band);
  TIFFGetField(band_red, TIFFTAG_SAMPLEFORMAT, &sample_bands);

  TIFFClose(band_red);

  // ================= COMPUTE REFLECTANCE & RADIANCE

  double radiance_matriz[8][height_band][width_band];
  double reflectance_matriz[8][height_band][width_band];

  std::ofstream outputReflectance("../out_data_txt/reflectance.txt"); 
  std::streambuf* coutReflectance = std::cout.rdbuf();
  std::cout.rdbuf(outputReflectance.rdbuf());

  for (int i = 1; i < 8; i++) {
    string path_tiff_band = argv[i];
    TIFF *curr_band = TIFFOpen(path_tiff_band.c_str(), "rm");

    tdata_t curr_band_line_buff = _TIFFmalloc(TIFFScanlineSize(curr_band));
    unsigned short curr_band_byte_size = TIFFScanlineSize(curr_band) / width_band;
    PixelReader curr_band_pixel_read = PixelReader(sample_bands, curr_band_byte_size, curr_band_line_buff);

    std::cout << "reflectance - band " << i << std::endl;
    for(int line = 0; line < height_band; line++){
        TIFFReadScanline(curr_band, curr_band_line_buff, line);

        for(int col = 0; col < width_band; col++){
            double band_pixel = curr_band_pixel_read.read_pixel(col);

            double radiance_pixel = band_pixel * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE];
            double reflectance_pixel = (PI * radiance_pixel * mtl.distance_earth_sun * mtl.distance_earth_sun) / (sensor.parameters[i][sensor.ESUN] * cos_theta);

            radiance_matriz[i][line][col] = radiance_pixel;
            reflectance_matriz[i][line][col] = reflectance_pixel;

            std::cout << reflectance_pixel << " ";
        }
        std::cout << std::endl;
    }

    _TIFFfree(curr_band_line_buff);
  }


  // ================= COMPUTE ALBEDO

  double albedo_matriz[height_band][width_band];

  std::ofstream outputAlbedo("../out_data_txt/albedo.txt"); 
  std::streambuf* coutAlbedo = std::cout.rdbuf();
  std::cout.rdbuf(outputAlbedo.rdbuf());

  std::cout << "albedo" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    for (int col = 0; col < width_band; col++) {

      double alb = 0;
      for (int i = 0; i < 7; i++) {
        alb = reflectance_matriz[1][line][col] * sensor.parameters[1][sensor.WB] +
            reflectance_matriz[2][line][col] * sensor.parameters[2][sensor.WB] +
            reflectance_matriz[3][line][col] * sensor.parameters[3][sensor.WB] +
            reflectance_matriz[4][line][col] * sensor.parameters[4][sensor.WB] +
            reflectance_matriz[5][line][col] * sensor.parameters[5][sensor.WB] +
            reflectance_matriz[7][line][col] * sensor.parameters[7][sensor.WB];
      }
      albedo_matriz[line][col] = alb;

      std::cout << alb << " ";
    }
   std::cout << std::endl;
  }


  // ================= COMPUTE NDVI

  double ndvi_matriz[height_band][width_band];

  std::ofstream outputNDVI("../out_data_txt/ndvi.txt"); 
  std::streambuf* coutNDVI = std::cout.rdbuf();
  std::cout.rdbuf(outputNDVI.rdbuf());

  std::cout << "NDVI" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    for (int col = 0; col < width_band; col++) {
      double ndvi = (reflectance_matriz[4][line][col] - reflectance_matriz[3][line][col]) / (reflectance_matriz[4][line][col] + reflectance_matriz[3][line][col]);
      ndvi_matriz[line][col] = ndvi;
      std::cout << ndvi << " ";
    }
   std::cout << std::endl;
  }


  // ================= COMPUTE LAI

  double lai_matriz[height_band][width_band];

  std::ofstream outputLAI("../out_data_txt/lai.txt"); 
  std::streambuf* coutLAI = std::cout.rdbuf();
  std::cout.rdbuf(outputLAI.rdbuf());

  std::cout << "LAI" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    double savi_line[width_band];
    double L = 0.05;

    for (int col = 0; col < width_band; col++){
      double savi = ((1 + L) * (reflectance_matriz[4][line][col] - reflectance_matriz[3][line][col])) / (L + (reflectance_matriz[4][line][col] + reflectance_matriz[3][line][col]));
      
      double lai;
      if (!isnan(savi) && definitelyGreaterThan(savi, 0.687))
          lai = 6;
      else if (!isnan(savi) && definitelyLessThan(savi, 0.1))
          lai = 0;
      else
          lai = -log((0.69 - savi) / 0.59) / 0.91;

      lai_matriz[line][col] = lai;
      std::cout << lai << " ";
    }
    std::cout << std::endl;
  }


  // ================= COMPUTE SURFACE TEMP

  double surface_temperature_matriz[height_band][width_band];

  double k1, k2;
  switch(mtl.number_sensor){
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

  std::ofstream outputSurfaceTemp("../out_data_txt/surface_temperature.txt"); 
  std::streambuf* coutSurfaceTemp = std::cout.rdbuf();
  std::cout.rdbuf(outputSurfaceTemp.rdbuf());

  std::cout << "surface_temperature" << std::endl;
  for(int line = 0; line < height_band; line++){ 

    int radiance_number = (mtl.number_sensor == 5)? 6: 7;
    for(int col = 0; col < width_band; col++) {

        double enb_emissivity = 0.97 + 0.0033 * lai_matriz[line][col];
        if(definitelyLessThan(ndvi_matriz[line][col], 0) || definitelyGreaterThan(lai_matriz[line][col], 2.99))
            enb_emissivity = 0.98;

        double surface_temperature = k2 / (log( (enb_emissivity * k1 / radiance_matriz[radiance_number][line][col]) + 1));    
        surface_temperature_matriz[line][col] = surface_temperature;

      std::cout << surface_temperature << " ";
    }    
    std::cout << std::endl;
  }


  // ================= COMPUTE NET RADIATION

  double net_radiation_matriz[height_band][width_band];

  string path_tiff_band_tal = argv[INPUT_BAND_TAL_INDEX];
  TIFF *band_tal = TIFFOpen(path_tiff_band_tal.c_str(), "rm");

  uint32 height_tal, width_tal;
  uint16 sample_tal_band;

  TIFFGetField(band_tal, TIFFTAG_IMAGELENGTH, &height_tal);
  TIFFGetField(band_tal, TIFFTAG_IMAGEWIDTH, &width_tal);
  TIFFGetField(band_tal, TIFFTAG_SAMPLEFORMAT, &sample_tal_band);

  tdata_t tal_band_line_buff = _TIFFmalloc(TIFFScanlineSize(band_tal));
  unsigned short tall_band_byte_size = TIFFScanlineSize(band_tal) / width_band;
  PixelReader tal_band_pixel_read = PixelReader(sample_tal_band, tall_band_byte_size, tal_band_line_buff);

  double temperature_kelvin = station.temperature_image + 273.15;
  double temperature_kelvin_pow_4 = temperature_kelvin * temperature_kelvin * temperature_kelvin * temperature_kelvin;

  std::ofstream outputRN("../out_data_txt/RN.txt"); 
  std::streambuf* coutRN = std::cout.rdbuf();
  std::cout.rdbuf(outputRN.rdbuf());
  std::cout << "RN" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    TIFFReadScanline(band_tal, tal_band_line_buff, line);

    for(int col = 0; col < width_band; col++){
      double tal_pixel = tal_band_pixel_read.read_pixel(col);
      double ea_emissivity = 0.85 * pow((-1 * log(tal_pixel)), 0.09); //TODO: Olhar esse cara ~ EA_terra

      double eo_emissivity = 0.95 + 0.01 * lai_matriz[line][col];
      if (definitelyLessThan(ndvi_matriz[line][col], 0) || definitelyGreaterThan(lai_matriz[line][col], 2.99))
        double eo_emissivity = 0.98;

      double temperature_pixel = surface_temperature_matriz[line][col];
      double surface_temperature_pow_4 = temperature_pixel * temperature_pixel * temperature_pixel * temperature_pixel;


      double short_wave_radiation = (1367 * cos_theta * tal_pixel) / (mtl.distance_earth_sun * mtl.distance_earth_sun);
      double large_wave_radiation_atmosphere = ea_emissivity * 5.67 * 1e-8 * temperature_kelvin_pow_4; // NaN
      double large_wave_radiation_surface = eo_emissivity * 5.67 * 1e-8 * surface_temperature_pow_4;


      double rn = short_wave_radiation - (short_wave_radiation * albedo_matriz[line][col]) +
                  large_wave_radiation_atmosphere - large_wave_radiation_surface -
                  (1 - eo_emissivity) * large_wave_radiation_atmosphere;

      if (definitelyLessThan(rn, 0)) 
        rn = 0;

      net_radiation_matriz[line][col] = rn;
      std::cout << rn << " ";
    }
   std::cout << std::endl;
  }

  TIFFClose(band_tal);


  // ================= COMPUTE SOIL HEAT (G)


  return 0;
}