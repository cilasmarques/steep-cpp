#include <fstream>
#include <iostream>

#include "landsat.h"
#include "constants.h"
#include "parameters.h"
#include "reader.h"

int main(int argc, char *argv[])
{
  int INPUT_BAND_BLUE_INDEX    = 1;
  int INPUT_BAND_GREEN_INDEX   = 2;
  int INPUT_BAND_RED_INDEX     = 3;
  int INPUT_BAND_NIR_INDEX     = 4;
  int INPUT_BAND_SWIR1_INDEX   = 5;
  int INPUT_BAND_TERMAL_INDEX  = 6;
  int INPUT_BAND_SWIR2_INDEX   = 7;
  int INPUT_BAND_TAL_INDEX     = 8;
  int INPUT_MTL_DATA_INDEX     = 9;
  int INPUT_STATION_DATA_INDEX = 10;
  int INPUT_LAND_COVER_INDEX   = 11;
  int OUTPUT_FOLDER            = 12;

  // load meta data
  string path_meta_file = argv[INPUT_MTL_DATA_INDEX];
  MTL mtl = MTL(path_meta_file);

  // load station data
  string station_data_path = argv[INPUT_STATION_DATA_INDEX];
  Station station = Station(station_data_path, mtl.image_hour);

  // load sensor data
  Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

  // load bands path
  string bands_paths[8];
  for (int i = 1; i < 8; i++) {
    bands_paths[i] = argv[i];
  }

  // load tal path
  string tal_path = argv[INPUT_BAND_TAL_INDEX];

  // load land cover path
  string land_cover_path = (argc >= 13) ? argv[12] : "";

  // load selected method 
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

  // load threads number
  int threads_num = 8;
  if(argc >= 16){
    string threads_flag = argv[15];
    if(threads_flag.substr(0,9) == "-threads=")
      threads_num = atof(threads_flag.substr(9, threads_flag.size()).c_str());
  }

  // load blocks number
  int blocks_num = 6504;
  if(argc >= 17){
    string blocks_flag = argv[16];
    if(blocks_flag.substr(0,8) == "-blocks=")
      blocks_num = atof(blocks_flag.substr(8, blocks_flag.size()).c_str());
  }

  //Timing
  using namespace std::chrono;
  system_clock::time_point begin, end;
  int64_t initial_time, final_time, general_time;
  std::cout << "PHASE,TIMESTAMP,START_TIME,END_TIME" << std::endl;

  initial_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  begin = system_clock::now();

  Landsat landsat = Landsat(method, bands_paths, tal_path, land_cover_path, threads_num, blocks_num);
  landsat.process_products(mtl, sensor, station);

  end = system_clock::now();
  final_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
  general_time = duration_cast<milliseconds>(end.time_since_epoch() - begin.time_since_epoch()).count();

  std::cout << "TOTAL," << general_time << "," << initial_time << "," << final_time << std::endl;

  return 0;
}
