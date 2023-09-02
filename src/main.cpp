#include <fstream>
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"
#include "STEEP.h"
#include "debug.h"
#include "products.h"

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


  std::ofstream output("./output/products.txt"); 
  std::streambuf* coutOutput = std::cout.rdbuf();
  std::cout.rdbuf(output.rdbuf());

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

  for (int i = 1; i < 8; i++) {
    string path_tiff_band = argv[i];
    TIFF *curr_band = TIFFOpen(path_tiff_band.c_str(), "rm");

    tdata_t curr_band_line_buff = _TIFFmalloc(TIFFScanlineSize(curr_band));
    unsigned short curr_band_byte_size = TIFFScanlineSize(curr_band) / width_band;
    PixelReader curr_band_pixel_read = PixelReader(sample_bands, curr_band_byte_size, curr_band_line_buff);

    std::cout << " ==== reflectance - band " << i << std::endl;
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
  int final_tif_calc = mtl.number_sensor == 8 ? 6 : 7;

  string tal_path = argv[INPUT_BAND_TAL_INDEX];
  TIFF *tal = TIFFOpen(tal_path.c_str(), "rm");

  tdata_t tal_line_buff = _TIFFmalloc(TIFFScanlineSize(tal));
  unsigned short curr_tal_line_size = TIFFScanlineSize(tal) / width_band;
  PixelReader tal_reader = PixelReader(sample_bands, curr_tal_line_size, tal_line_buff);

  std::cout << " ==== albedo" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    TIFFReadScanline(tal, tal_line_buff, line);

    for (int col = 0; col < width_band; col++) {

      double alb = 0;
      for (int i = 0; i < 7; i++) {
        alb = reflectance_matriz[1][line][col] * sensor.parameters[1][sensor.WB] +
            reflectance_matriz[2][line][col] * sensor.parameters[2][sensor.WB] +
            reflectance_matriz[3][line][col] * sensor.parameters[3][sensor.WB] +
            reflectance_matriz[4][line][col] * sensor.parameters[4][sensor.WB] +
            reflectance_matriz[5][line][col] * sensor.parameters[5][sensor.WB] +
            reflectance_matriz[final_tif_calc][line][col] * sensor.parameters[final_tif_calc][sensor.WB];
      }

      alb = (alb - 0.03) / (tal_reader.read_pixel(col) * tal_reader.read_pixel(col));

      albedo_matriz[line][col] = alb;

      std::cout << alb << " ";
    }

   std::cout << std::endl;
  }
  _TIFFfree(tal_line_buff);
  TIFFClose(tal);

  // ================= COMPUTE NDVI

  double ndvi_matriz[height_band][width_band];

  std::cout << " ==== NDVI" << std::endl;
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

  std::cout << " ==== LAI" << std::endl;
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


  // ================= COMPUTE PAI

  double pai_matriz[height_band][width_band];

  std::cout << " ==== PAI" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    for (int col = 0; col < width_band; col++){
      double pai_value = 10.1 * (reflectance_matriz[4][line][col] - sqrt(reflectance_matriz[3][line][col])) + 3.1;

      if (pai_value < 0) {
        pai_value = 0;
      }

      std::cout << pai_value << " ";
      pai_matriz[line][col] = pai_value;
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

  std::cout << " ==== surface_temperature" << std::endl;
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

  std::cout << " ==== RN" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    TIFFReadScanline(band_tal, tal_band_line_buff, line);

    for(int col = 0; col < width_band; col++){
      double tal_pixel = tal_band_pixel_read.read_pixel(col);
      double ea_emissivity = 0.85 * pow((-1 * log(tal_pixel)), 0.09);
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

  double G_matriz[height_band][width_band];

  std::cout << " ==== G - soil heat flux" << std::endl;
  for(int line = 0; line < height_band; line++){ 
    for (int col = 0; col < width_band; col++)
    {
      double G;
      if (essentiallyEqual(ndvi_matriz[line][col], 0) || definitelyGreaterThan(ndvi_matriz[line][col], 0))
      {
        double ndvi_pixel_pow_4 = ndvi_matriz[line][col] * ndvi_matriz[line][col] * ndvi_matriz[line][col] * ndvi_matriz[line][col];
        G = (surface_temperature_matriz[line][col] - 273.15) * (0.0038 + 0.0074 * albedo_matriz[line][col]) *
                              (1 - 0.98 * ndvi_pixel_pow_4) * net_radiation_matriz[line][col];
      }
      else
        G = 0.5 * net_radiation_matriz[line][col];

      if (definitelyLessThan(G, 0))
        G = 0;
 
      G_matriz[line][col] = G;
      std::cout << G << " ";
    }
    std::cout << std::endl;
  }

  // ============== FINALIZA A PRIMEIRA ETAPA

  std::vector<std::vector<double>> pai_vector;
  std::vector<std::vector<double>> ndvi_vector;
  std::vector<std::vector<double>> surface_temperature_vector;
  std::vector<std::vector<double>> albedo_vector;
  std::vector<std::vector<double>> net_radiation_vector;
  std::vector<std::vector<double>> soil_heat_vector;

  for (int i = 0; i < height_band; i++) {
    std::vector<double> pai_vector_interno;
    std::vector<double> ndvi_vector_interno;
    std::vector<double> surface_temperature_vector_interno;
    std::vector<double> albedo_vector_interno;
    std::vector<double> net_radiation_vector_interno;
    std::vector<double> soil_heat_vector_interno;
  
    for (int j = 0; j < width_band; j++) {
      pai_vector_interno.push_back(pai_matriz[i][j]);
      ndvi_vector_interno.push_back(ndvi_matriz[i][j]);
      surface_temperature_vector_interno.push_back(surface_temperature_matriz[i][j]);
      albedo_vector_interno.push_back(albedo_matriz[i][j]);
      net_radiation_vector_interno.push_back(net_radiation_matriz[i][j]);
      soil_heat_vector_interno.push_back(G_matriz[i][j]);
    }

    pai_vector.push_back(pai_vector_interno);
    ndvi_vector.push_back(ndvi_vector_interno);
    surface_temperature_vector.push_back(surface_temperature_vector_interno);
    albedo_vector.push_back(albedo_vector_interno);
    net_radiation_vector.push_back(net_radiation_vector_interno);
    soil_heat_vector.push_back(soil_heat_vector_interno);
  }


  // ============== SELEÇÃO DE PIXEL

  Candidate hot_pixel, cold_pixel;
  hot_pixel = getHotPixelSTEPP(ndvi_vector, surface_temperature_vector, albedo_vector, net_radiation_vector, soil_heat_vector, height_band, width_band);
  cold_pixel = getColdPixelSTEPP(ndvi_vector, surface_temperature_vector, albedo_vector, net_radiation_vector, soil_heat_vector, height_band, width_band);


  // ============== CALCULO DO CICLO RAH INICIAL

  // talvez tenha algo de errado com esses caras
  double ustar_station = (VON_KARMAN * station.v6) / (log(station.WIND_SPEED / station.SURFACE_ROUGHNESS));
  double u10 = (ustar_station / VON_KARMAN) * log(10 / station.SURFACE_ROUGHNESS); // TODO: Talvez esse cara esteja errado

  double ndvi_min = 1.0;
  double ndvi_max = 0.0;
  for (int line = 0; line < height_band; line++)
  {
    std::vector<double> ndvi_line = ndvi_vector[line];
    for (int col = 0; col < width_band; col++)
    {
      if (ndvi_line[col] < ndvi_min)
        ndvi_min = ndvi_line[col];
      if (ndvi_line[col] > ndvi_max)
        ndvi_max = ndvi_line[col];
    }
  }

  double d0_matriz[height_band][35];
  double zom_matriz[height_band][35];
  double ustar_matriz[height_band][35];
  double kb1_matriz[height_band][35];
  double aerodynamic_resistance_matriz[height_band][35];

  for (int line = 0; line < height_band; line++)
  {
    d0_fuction(pai_vector[line], width_band, d0_matriz[line]);
    zom_fuction(station.A_ZOM, station.B_ZOM, ndvi_matriz[line], width_band, zom_matriz[line], d0_matriz[line], pai_vector[line]);
    ustar_fuction(u10, zom_matriz[line], width_band, ustar_matriz[line], d0_matriz[line]);
    kb_function(ustar_matriz[line], zom_matriz[line], pai_vector[line], ndvi_matriz[line], ndvi_max, ndvi_min, width_band, kb1_matriz[line]);
    aerodynamic_resistance_fuction(ustar_matriz[line], width_band, aerodynamic_resistance_matriz[line], zom_matriz[line], d0_matriz[line], kb1_matriz[line]);
  }

  std::vector<std::vector<double>> d0_vector = arrayToVector(d0_matriz, height_band, width_band);
  std::vector<std::vector<double>> zom_vector = arrayToVector(zom_matriz, height_band, width_band);
  std::vector<std::vector<double>> ustar_vector = arrayToVector(ustar_matriz, height_band, width_band);
  std::vector<std::vector<double>> kb1_vector = arrayToVector(kb1_matriz, height_band, width_band);
  std::vector<std::vector<double>> aerodynamic_resistance_vector = arrayToVector(aerodynamic_resistance_matriz, height_band, width_band);

  std::cout << " ==== d0" << std::endl;
  printVector2x2(d0_vector);

  std::cout << " ==== zom" << std::endl;
  printVector2x2(zom_vector);

  std::cout << " ==== ustar" << std::endl;
  printVector2x2(ustar_vector);

  std::cout << " ==== kb1" << std::endl;
  printVector2x2(kb1_vector);

  std::cout << " ==== rah" << std::endl;
  printVector2x2(aerodynamic_resistance_vector);

  // ============== CALCULO DO CICLO RAH FINAL

  // Auxiliaries arrays calculation
  double L[width_band];
  double y_01_line[width_band], y_2_line[width_band], x_200_line[width_band];
  double psi_01_line[width_band], psi_2_line[width_band], psi_200_line[width_band];

  double sensible_heat_flux_matriz[height_band][35];
  std::vector<double> ustar_write_line, aerodynamic_resistance_write_line;
  std::vector<std::vector<double>> ustar_previous, aerodynamic_resistance_previous;

  // pega a resistencia aerodinamica do PQ e adiciona no objeto hot_pixel
  double hot_pixel_aerodynamic = aerodynamic_resistance_vector[hot_pixel.line][hot_pixel.col];
  hot_pixel.aerodynamic_resistance.push_back(hot_pixel_aerodynamic);

  // pega a resistencia aerodinamica do PF e adiciona no objeto cold_pixel
  double cold_pixel_aerodynamic = aerodynamic_resistance_vector[cold_pixel.line][cold_pixel.col];
  cold_pixel.aerodynamic_resistance.push_back(cold_pixel_aerodynamic);

  // calcula o FC para o pixel quente e para o pixel frio
  double fc_hot = 1 - pow((ndvi_vector[hot_pixel.line][hot_pixel.col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);
  double fc_cold = 1 - pow((ndvi_vector[cold_pixel.line][cold_pixel.col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);

  double rah_aux_hot;
  double rah_aux_cold;

  double H_pf_terra;
  double H_pq_terra;

  double rah_ini_pq_terra;
  double rah_ini_pf_terra;

  for (int i = 0; i < 2; i++)
  {
    ustar_previous = ustar_vector;                                    //aux
    aerodynamic_resistance_previous = aerodynamic_resistance_vector;  //aux

    rah_ini_pq_terra = hot_pixel.aerodynamic_resistance[i];
    rah_ini_pf_terra = cold_pixel.aerodynamic_resistance[i];

    double LEc_terra = 0.55 * fc_hot * (hot_pixel.net_radiation - hot_pixel.soil_heat_flux) * 0.78;
    double LEc_terra_pf = 1.75 * fc_cold * (cold_pixel.net_radiation - cold_pixel.soil_heat_flux) * 0.78;

    H_pf_terra = cold_pixel.net_radiation - cold_pixel.soil_heat_flux - LEc_terra_pf;
    double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);

    H_pq_terra = hot_pixel.net_radiation - hot_pixel.soil_heat_flux - LEc_terra;
    double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);

    double x = aerodynamic_resistance_previous[hot_pixel.line][hot_pixel.col];

    double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
    double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));        

    std::cout << " ==== H ini " << i << std::endl;
    for (int line = 0; line < height_band; line++)
    {
      for (int col = 0; col < width_band; col++)
      {
        double DISP = d0_matriz[line][col];
        double dT_ini_terra = (a + b * (surface_temperature_vector[line][col] - 273.15));
        
        // H_ini_terra
        sensible_heat_flux_matriz[line][col] = RHO * SPECIFIC_HEAT_AIR * (dT_ini_terra) / aerodynamic_resistance_previous[line][col]; 

        // L_MB_terra
        double ustar_pow_3 = ustar_previous[line][col] * ustar_previous[line][col] * ustar_previous[line][col];
        L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_vector[line][col]) / (VON_KARMAN * GRAVITY * sensible_heat_flux_matriz[line][col]));

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

        // u_ast_corr_terra
        double ust = (VON_KARMAN * ustar_previous[line][col]) / (log((10 - DISP) / zom_matriz[line][col]) - psi_200_line[col]);
        ustar_write_line.push_back(ust);

        double zoh_terra = zom_matriz[line][col] / pow(exp(1.0), (kb1_matriz[line][col]));
        double temp_rah1_corr_terra = (ust * VON_KARMAN);
        double temp_rah2_corr_terra = log((10 - DISP) / zom_matriz[line][col]) - psi_2_line[col];
        double temp_rah3_corr_terra = temp_rah1_corr_terra * log(zom_matriz[line][col] / zoh_terra);        
        double rah = (temp_rah1_corr_terra * temp_rah2_corr_terra) + temp_rah3_corr_terra;

        aerodynamic_resistance_write_line.push_back(rah);

        if (line == hot_pixel.line && col == hot_pixel.col)
        {
          rah_aux_hot = aerodynamic_resistance_write_line[col];
          hot_pixel.aerodynamic_resistance.push_back(rah_aux_hot);
        }

        if (line == cold_pixel.line && col == cold_pixel.col)
        {
          rah_aux_cold = aerodynamic_resistance_write_line[col];
          cold_pixel.aerodynamic_resistance.push_back(rah_aux_cold);
        }

        std::cout << sensible_heat_flux_matriz[line][col] << " ";
      }
      std::cout << std::endl;

      ustar_vector[line] = ustar_write_line;
      aerodynamic_resistance_vector[line] = aerodynamic_resistance_write_line;

      ustar_write_line.clear();
      aerodynamic_resistance_write_line.clear();
    }
  }

  std::cout << " ==== rah final" << std::endl;
  printVector2x2(aerodynamic_resistance_vector);


  // ============== CALCULO DA EVOAPOTRANSPITAÇÃO

  double dt_pq_terra = H_pq_terra * rah_ini_pq_terra / (RHO * SPECIFIC_HEAT_AIR);
  double dt_pf_terra = H_pf_terra * rah_ini_pf_terra / (RHO * SPECIFIC_HEAT_AIR);

  double b = (dt_pq_terra - dt_pf_terra) / (hot_pixel.temperature - cold_pixel.temperature);
  double a = dt_pf_terra - (b * (cold_pixel.temperature - 273.15));

  std::cout << " ==== H final" << std::endl;
  // Sensible heat flux H
  for (int line = 0; line < height_band; line++) {
    for (int col = 0; col < width_band; col++) {
      sensible_heat_flux_matriz[line][col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_vector[line][col] - 273.15)) / aerodynamic_resistance_vector[line][col];

      if (!isnan(sensible_heat_flux_matriz[line][col]) && sensible_heat_flux_matriz[line][col] > (net_radiation_vector[line][col] - soil_heat_vector[line][col])) {
        sensible_heat_flux_matriz[line][col] = net_radiation_vector[line][col] - soil_heat_vector[line][col];
      }

      std::cout << sensible_heat_flux_matriz[line][col] << " ";
    }
    std::cout << std::endl;
  }


  double latent_heat_flux_matriz[height_band][35];
  double net_radiation_24h_matriz[height_band][35];
  double evapotranspiration_fraction_matriz[height_band][35];
  double sensible_heat_flux_24h_matriz[height_band][35];
  double latent_heat_flux_24h_matriz[height_band][35];
  double evapotranspiration_24h_matriz[height_band][35];
  double evapotranspiration_ulisses_matriz[height_band][35];

  //Upscalling temporal
  double dr = (1 / mtl.distance_earth_sun) * (1 / mtl.distance_earth_sun);
  double sigma = 0.409*sin(((2*PI/365)*mtl.julian_day)-1.39);
  double phi = (PI/180) * station.latitude;
  double omegas = acos(-tan(phi)*tan(sigma));
  double Ra24h = (((24*60/PI)*GSC*dr)*(omegas*sin(phi)*
                              sin(sigma)+cos(phi)*cos(sigma)*sin(omegas)))*(1000000/86400.0);

  //Short wave radiation incident in 24 hours (Rs24h)
  double Rs24h = station.INTERNALIZATION_FACTOR * sqrt(station.v7_max - station.v7_min) * Ra24h;
  
  // Latent heat flux LE & Compute evapotranspiration
  for(int line = 0; line < height_band; line++){
      latent_heat_flux_function(net_radiation_matriz[line], G_matriz[line], sensible_heat_flux_matriz[line], width_band, latent_heat_flux_matriz[line]);
      net_radiation_24h_function(albedo_matriz[line], Ra24h, Rs24h, width_band, net_radiation_24h_matriz[line]);
      evapotranspiration_fraction_fuction(latent_heat_flux_matriz[line], net_radiation_matriz[line], G_matriz[line], width_band, evapotranspiration_fraction_matriz[line]);
      sensible_heat_flux_24h_fuction(evapotranspiration_fraction_matriz[line], net_radiation_24h_matriz[line], width_band, sensible_heat_flux_24h_matriz[line]);
      latent_heat_flux_24h_function(evapotranspiration_fraction_matriz[line], net_radiation_24h_matriz[line], width_band, latent_heat_flux_24h_matriz[line]);
      evapotranspiration_24h_function(latent_heat_flux_24h_matriz[line], station, width_band, evapotranspiration_24h_matriz[line]);
      evapotranspiration_ulisses_function(net_radiation_24h_matriz[line], evapotranspiration_fraction_matriz[line], width_band, evapotranspiration_ulisses_matriz[line]);
  }

  std::vector<std::vector<double>> latent_heat_flux_vector = arrayToVector(latent_heat_flux_matriz, height_band, width_band);
  std::vector<std::vector<double>> net_radiation_24h_vector = arrayToVector(net_radiation_24h_matriz, height_band, width_band);
  std::vector<std::vector<double>> evapotranspiration_fraction_vector = arrayToVector(evapotranspiration_fraction_matriz, height_band, width_band);
  std::vector<std::vector<double>> sensible_heat_flux_24h_vector = arrayToVector(sensible_heat_flux_24h_matriz, height_band, width_band);
  std::vector<std::vector<double>> latent_heat_flux_24h_vector = arrayToVector(latent_heat_flux_24h_matriz, height_band, width_band);
  std::vector<std::vector<double>> evapotranspiration_24h_vector = arrayToVector(evapotranspiration_24h_matriz, height_band, width_band);
  std::vector<std::vector<double>> evapotranspiration_ulisses_vector = arrayToVector(evapotranspiration_ulisses_matriz, height_band, width_band);

  std::cout << " ==== LE" << std::endl;
  printVector2x2(latent_heat_flux_vector);

  std::cout << " ==== Rn24" << std::endl;
  printVector2x2(net_radiation_24h_vector);

  std::cout << " ==== evapotranspiration_fraction" << std::endl;
  printVector2x2(evapotranspiration_fraction_vector);

  std::cout << " ==== H24" << std::endl;
  printVector2x2(sensible_heat_flux_24h_vector);

  std::cout << " ==== LE24" << std::endl;
  printVector2x2(latent_heat_flux_24h_vector);

  std::cout << " ==== ET24" << std::endl;
  printVector2x2(evapotranspiration_24h_vector);

  std::cout << " ==== ET" << std::endl;
  printVector2x2(evapotranspiration_ulisses_vector);

  return 0;
}