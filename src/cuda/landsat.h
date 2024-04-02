#pragma once

#include "products.h"

#include "debug.h"
#include "utils.h"
#include "reader.h"
#include "constants.h"
#include "candidate.h"
#include "constants.h"
#include "endmembers.h"
#include "parameters.h"

/**
 * @brief  Struct to manage the products calculation.
 */
struct Landsat
{
  int method;
  int threads_num;
  int blocks_num;

  string tal_path;
  string land_cover_path;
  string bands_paths[8];

  Landsat(int method, string bands_paths[], string tal_path, string land_cover_path, int threads_num, int blocks_num);

  void process_products(MTL mtl, Sensor sensor, Station station);
};