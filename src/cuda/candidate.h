#pragma once

#include "utils.h"

/**
 * @brief  Struct representing a hot or cold pixel candidate.
 */
struct Candidate
{
  int line, col;
  int negative_neighbour;
  double coefficient_variation;
  double ndvi, temperature, ustar;
  double net_radiation, soil_heat_flux, ho, zom;
  vector<double> aerodynamic_resistance;

  /**
   * @brief  Empty constructor, all attributes are initialized with 0.
   */
  Candidate();

  /**
   * @brief  Copy constructor.
   * @param  c: Candidate to be copied.
  */
  Candidate(const Candidate &c);

  /**
   * @brief  Constructor with initialization values to attributes.
   * @param  ndvi: Pixel's NDVI.
   * @param  temperature: Pixel's surface temperature.
   * @param  net_radiation: Pixel's net radiation.
   * @param  soil_heat_flux: Pixel's soil heat flux.
   * @param  ho: Pixel's ho.
   * @param  line: Pixel's line on TIFF.
   * @param  col: Pixel's column on TIFF.
   */
  Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho, int line, int col);

  /**
   * @brief  Calculates a initial value for Pixel's aerodynamic resistance. Adding this value to attribute aerodynamic resistance.
   * @param  u200: Wind speed at 200 m.
   * @param  A_ZOM: Coefficient A.
   * @param  B_ZOM: Coefficient B.
   * @param  VON_KARMAN: Karman's constant.
   */
  void setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN);

  /**
   * @brief  Prints the data contained at the struct.
   */
  void toString();
};

/**
 * @brief  Compares two Candidates based upon their position.
 * @param  a: First candidate.
*/
bool equals(Candidate a, Candidate b);

/**
 * @brief  Compares two Candidates based upon their surface temperature.
 * @param  a: First candidate.
 * @param  b: Second candidate.
 * @retval TRUE if second candidate is greater than first one, and FALSE otherwise.
 */
bool compare_candidate_temperature(Candidate a, Candidate b);

/**
 * @brief  Compares two Candidates based upon their NDVI.
 * @param  a: First candidate.
 * @param  b: Second candidate.
 * @retval TRUE if second candidate is greater than first one, and FALSE otherwise.
 */
bool compare_candidate_ndvi(Candidate a, Candidate b);

/**
 * @brief  Compares two Candidates based upon their HO.
 * @param  a: First candidate.
 * @param  b: Second candidate.
 * @retval TRUE if second candidate is greater than first one, and FALSE otherwise.
 */
bool compare_candidate_ho(Candidate a, Candidate b);