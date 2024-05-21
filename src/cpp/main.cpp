#include <fstream>
#include <iostream>

#include "utils.h"
#include "reader.h"
#include "landsat.h"
#include "constants.h"
#include "parameters.h"

/** 
 * @brief Main function
 * This function is responsible for reading the input parameters and calling the Landsat class to process the products.
 * 
 * @param argc Number of input parameters
 * @param argv Input parameters
 *              - INPUT_BAND_BLUE_INDEX    = 1;
 *              - INPUT_BAND_GREEN_INDEX   = 2;
 *              - INPUT_BAND_RED_INDEX     = 3;
 *              - INPUT_BAND_NIR_INDEX     = 4;
 *              - INPUT_BAND_SWIR1_INDEX   = 5;
 *              - INPUT_BAND_TERMAL_INDEX  = 6;
 *              - INPUT_BAND_SWIR2_INDEX   = 7;
 *              - INPUT_BAND_TAL_INDEX     = 8;
 *              - INPUT_MTL_DATA_INDEX     = 9;
 *              - INPUT_STATION_DATA_INDEX = 10;
 *              - INPUT_LAND_COVER_INDEX   = 11;
 *              - OUTPUT_FOLDER            = 12;
 * @return int
*/
int main(int argc, char *argv[])
{
  int INPUT_BAND_TAL_INDEX    = 8;
  int INPUT_MTL_DATA_INDEX     = 9;
  int INPUT_STATION_DATA_INDEX = 10;
  int INPUT_LAND_COVER_INDEX   = 11;
  int OUTPUT_FOLDER            = 12;
  int METHOD_INDEX             = 13;
  int THREADS_INDEX            = 14;
  int BLOCKS_INDEX             = 15;

  string path_meta_file = argv[INPUT_MTL_DATA_INDEX];
  string station_data_path = argv[INPUT_STATION_DATA_INDEX];
  string land_cover_path = argv[INPUT_LAND_COVER_INDEX];

  // load bands path
  string bands_paths[INPUT_BAND_TAL_INDEX + 1];
  for (int i = 1; i <= INPUT_BAND_TAL_INDEX; i++) {
    bands_paths[i] = argv[i];
  }

  // load selected method 
  int method = 0;
  if(argc >= 14){
    string flag = argv[METHOD_INDEX];
    if(flag.substr(0, 6) == "-meth=")
      method = flag[6] - '0';
  }

  // load threads number
  int threads_num = 1;
  if(argc >= 15){
    string threads_flag = argv[THREADS_INDEX];
    if(threads_flag.substr(0,9) == "-threads=")
      threads_num = atof(threads_flag.substr(9, threads_flag.size()).c_str());
  }

  // load output folder
  string output_folder = argv[OUTPUT_FOLDER];
  string output_time = output_folder + "/time.csv";  
  string output_metadata = output_folder + "/metadata.txt";
  string output_products = output_folder + "/products.txt";

  // =====  START + TIME OUTPUT =====
  MTL mtl = MTL(path_meta_file);
  Sensor sensor = Sensor(mtl.number_sensor, mtl.year);
  Station station = Station(station_data_path, mtl.image_hour);
  Landsat landsat = Landsat(bands_paths, land_cover_path, mtl);

  ofstream time_output;
  time_output.open(output_time);
  string general, initial, final;
  system_clock::time_point begin, end;
  int64_t initial_time, final_time, general_time;

  begin = system_clock::now();
  initial_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();

  time_output << "PHASE,TIMESTAMP,START_TIME,END_TIME" << std::endl;
  time_output << landsat.compute_Rn_G(sensor, station);
  time_output << landsat.select_endmembers(method);
  time_output << landsat.converge_rah_cycle(station, method, threads_num);
  time_output << landsat.compute_H_ET(station);

  end = system_clock::now();
  final_time = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
  general_time = duration_cast<nanoseconds>(end.time_since_epoch() - begin.time_since_epoch()).count();
  time_output << "TOTAL," << general_time << "," << initial_time << "," << final_time << std::endl;
  time_output.close();

  landsat.save_products(output_products);
  landsat.close();

  return 0;
}
