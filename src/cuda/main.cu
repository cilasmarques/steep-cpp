#include <fstream>
#include <iostream>

#include "landsat.h"
#include "constants.h"
#include "parameters.h"
#include "reader.h"

using namespace std;

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

  // load meta data
  string path_meta_file = argv[INPUT_MTL_DATA_INDEX];
  MTL mtl = MTL(path_meta_file);

  // load station data
  string station_data_path = argv[INPUT_STATION_DATA_INDEX];
  Station station = Station(station_data_path, mtl.image_hour);

  // load sensor data
  Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

  // load bands path
  string bands_paths[INPUT_BAND_TAL_INDEX];
  for (int i = 1; i < INPUT_BAND_TAL_INDEX; i++) {
    bands_paths[i] = argv[i];
  }

  // load tal path
  string tal_path = argv[INPUT_BAND_TAL_INDEX];

  // load land cover path
  string land_cover_path = argv[INPUT_LAND_COVER_INDEX];

  // load selected method 
  int method = 0;
  if(argc >= 14){
    string flag = argv[13];
    if(flag.substr(0, 6) == "-meth=")
      method = flag[6] - '0';
  }

  // load threads number
  int threads_num = 8;
  if(argc >= 15){
    string threads_flag = argv[14];
    if(threads_flag.substr(0,9) == "-threads=")
      threads_num = atof(threads_flag.substr(9, threads_flag.size()).c_str());
  }

  // load blocks number
  int blocks_num = 6504;
  if(argc >= 16){
    string blocks_flag = argv[15];
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

  // =====  END + OUTPUTS =====

  string outPath = std::string(argv[OUTPUT_FOLDER]) + "/metadata.txt";
  std::ofstream outputThreads(outPath);
  std::streambuf* coutThreads = std::cout.rdbuf();
  std::cout.rdbuf(outputThreads.rdbuf());
  std::cout << outPath << std::endl;
  std::cout << "informed threads: " << threads_num << std::endl;
  std::cout << "informed blocks: " << blocks_num << std::endl;

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device);
      
      std::cout << "Device " << device << " is a " << deviceProp.name << std::endl;
      std::cout << "Device " << device << " has " << deviceProp.multiProcessorCount << " SMs" << std::endl;
      std::cout << "Device " << device << " has " << deviceProp.persistingL2CacheMaxSize  << " bytes of L2 cache" << std::endl;
      std::cout << "Device " << device << " has " << deviceProp.concurrentKernels << " concurrent kernels" << std::endl;
      std::cout << "Device " << device << " has " << deviceProp.maxBlocksPerMultiProcessor << " max blocks per SM" << std::endl;
      std::cout << "Device " << device << " has " << deviceProp.maxThreadsPerMultiProcessor << " max threads per SM" << std::endl;
      std::cout << "Device " << device << " has " << deviceProp.maxGridSize << " max grid size" << std::endl;
  }

  return 0;
}
