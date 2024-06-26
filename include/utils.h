#pragma once

#include "constants.h"

/**
 * @brief  Determines if a and b are approximately equals based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if they are approximately equals, and FALSE otherwise.
 */
bool approximatelyEqual(double a, double b);

/**
 * @brief  Determines if a and b are essentially equals based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if they are essentially equals, and FALSE otherwise.
 */
bool essentiallyEqual(double a, double b);

/**
 * @brief  Determines if a is definitely greater than b based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if a is definitely greater than b, and FALSE otherwise.
 */
bool definitelyGreaterThan(double a, double b);

/**
 * @brief  Determines if a is definitely less than b based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if a is definitely less than b, and FALSE otherwise.
 */
bool definitelyLessThan(double a, double b);

/**
 * @brief  Prints a 2x2 matrix.
 * 
 * @param  matrix: Matrix to be printed.
 */
void printVector2x2(const vector<vector<double>> &matrix);
