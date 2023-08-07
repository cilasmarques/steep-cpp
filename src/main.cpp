#include <fstream>
#include <iostream>

#include "landsat.h"
#include "parameters.h"

// #include "candidate.h"
// #include "STEEP.h"
// #include "debug.h"
// #include "products.h"

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
  const int INPUT_LAND_COVER_INDEX = 11;
  const int OUTPUT_FOLDER = 12;
  // const int INPUT_ALG_METHOD_INDEX = 13;
  // const int INPUT_NAN_VALUE_INDEX = 14;


  // ================= LOAD DATA

  // load meta data
  string path_meta_file = argv[INPUT_MTL_DATA_INDEX];
  MTL mtl = MTL(path_meta_file);

  // load sensor data
  string station_data_path = argv[INPUT_STATION_DATA_INDEX];
  Station station = Station(station_data_path, mtl.image_hour);

  // load sensor data
  Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

  // load land cover data
  string landCoverPath = (argc >= 13) ? argv[12] : "";

  // load method data
  int method = 0;
  if(argc >= 14){
    string flag = argv[13];
    if(flag.substr(0, 6) == "-meth=")
      method = flag[6] - '0';
  }

  // load nan value
  double noData = NaN;
  if(argc >= 15){
    string noData_flag = argv[14];
    if(noData_flag.substr(0,5) == "-nan=")
      noData = atof(noData_flag.substr(5, noData_flag.size()).c_str());
  }

  // Open bands
  TIFF *bands_resampled[8];
  for (int i = 1; i < 8; i++)
  {
    string path_tiff_base = argv[i];
    bands_resampled[i] = TIFFOpen(path_tiff_base.c_str(), "rm");
    check_open_tiff(bands_resampled[i]);
  }

  string output_path = argv[OUTPUT_FOLDER];
  string tal_path = argv[INPUT_BAND_TAL_INDEX];

  Landsat landsat = Landsat(tal_path, output_path, method, noData, landCoverPath);

  landsat.process_partial_products(bands_resampled, mtl, station, sensor);

  return 0;
}
