#pragma once

#include "constants.h"
#include "utils.h"
#include "reader.h"
#include "parameters.h"
#include "candidate.h"

/**
 * @brief   Checks if the input value matches one of the agricultural land cover code.
 * @param   value: The input value.
 * @retval  TRUE if the values matches and FALSE otherwise.
 */
bool checkLandCode(int value);

/**
 * @brief   Tests the land cover homogeneity of a pixel.
 * @note    A the land cover of a pixel is homogeneous if every neighbour pixel inside a 7x7 window is also agricultural field.
 * @param   landCover: Land Cover TIFF.
 * @param   mask: Output binary TIFF, where pixels with 1 means that is a homogeneous pixel and 0 means the otherwise.
 */
void testLandCoverHomogeneity(TIFF* landCover, TIFF* mask);

/**
 * @brief   Tests the ndvi, surface_temperature and albedo homogeneity of a pixel.
 * @note    A pixel is homogeneous in these criteria if inside a 7x7 window the coefficient of variation of the albedo and ndvi is less or equal than 25%
 *          and the surface temperature has a standard deviation less or equal than 1.5 K.
 * @param   ndvi: NDVI TIFF.
 * @param   surface_temperature: TS TIFF.
 * @param   albedo: Albedo TIFF.
 * @param   maskLC: A binary TIFF conteining the data of the land cover homogeneity.
 * @param   output: A binary TIFF, where pixels with 1 means that is a homogeneous pixel in land cover, ndvi, surface temperature and albedo, and 0 means otherwise.
 */
void testHomogeneity(TIFF* ndvi, TIFF* surface_temperature, TIFF* albedo, TIFF* maskLC, TIFF* output);

/**
 * @brief   Removes from the binary TIFF gave as input, groups of pixel with value 1, that have less pixels than a specified value.
 * @param   input: A binary TIFF to be processed.
 * @param   output: A binary TIFF.
 * @param   groupSize: The inferior limit of the group of pixels size.
 */
void testMorphological(TIFF* input, TIFF* output, int groupSize);
