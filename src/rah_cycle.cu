#include "rah_cycle.cuh"

__global__ void correctionCycle(double *surfaceTemperatureLine, double *zomLine, double *ustarRLine, double *ustarWLine, double *rahRLine,
                                double *rahWLine, double *a, double *b, double *u200, int *size)
{

  // Identify position
  int pos = threadIdx.x + blockIdx.x * blockDim.x;

  printf("%d\n", size);

  while (pos < *size)
  {

    double sensibleHeatFlux = RHO * SPECIFIC_HEAT_AIR * (*a + *b * (surfaceTemperatureLine[pos] - 273.15)) / rahRLine[pos];

    double L = -1 * ((RHO * SPECIFIC_HEAT_AIR * pow(ustarRLine[pos], 3) * surfaceTemperatureLine[pos]) / (VON_KARMAN * GRAVITY * sensibleHeatFlux));

    double y01 = pow((1 - (16 * 0.1) / L), 0.25);
    double y2 = pow((1 - (16 * 2) / L), 0.25);
    double x200 = pow((1 - (16 * 200) / L), 0.25);

    double psi01, psi2, psi200;

    if (!isnan(L) && L > 0)
    {

      psi01 = -5 * (0.1 / L);
      psi2 = -5 * (2 / L);
      psi200 = -5 * (2 / L);
    }
    else
    {

      psi01 = 2 * log((1 + y01 * y01) / 2);

      psi2 = 2 * log((1 + y2 * y2) / 2);

      psi200 = 2 * log((1 + x200) / 2) + log((1 + x200 * x200) / 2) - 2 * atan(x200) + 0.5 * M_PI;
    }

    ustarWLine[pos] = (VON_KARMAN * *u200) / (log(200 / zomLine[pos]) - psi200);

    rahWLine[pos] = (log(2 / 0.1) - psi2 + psi01) / (ustarWLine[pos] * VON_KARMAN);

    pos += blockDim.x * gridDim.x;
  }
}

// TODO: Update that
// __global__ void correctionCycleSTEEP(double *surfaceTemperatureLine, double *zomLine, double *ustarRLine, double *ustarWLine, double *rahRLine,
//                                      double *rahWLine, double *a, double *b, double *u200, int *size)
// {

//   for (int col = 0; col < width_band; col++)
//   {
//     double DISP = d0_vector[line][col];
//     double dT_ini_terra = (a + b * (surface_temperature_vector[line][col] - 273.15));

//     // H_ini_terra
//     sensible_heat_flux_vector[line][col] = RHO * SPECIFIC_HEAT_AIR * (dT_ini_terra) / aerodynamic_resistance_previous[line][col];

//     // L_MB_terra
//     double ustar_pow_3 = ustar_previous[line][col] * ustar_previous[line][col] * ustar_previous[line][col];
//     L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_vector[line][col]) / (VON_KARMAN * GRAVITY * sensible_heat_flux_vector[line][col]));

//     y_01_line[col] = pow((1 - (16 * 0.1) / L[col]), 0.25);
//     y_2_line[col] = pow((1 - (16 * (10 - DISP)) / L[col]), 0.25);
//     x_200_line[col] = pow((1 - (16 * (10 - DISP)) / L[col]), 0.25);

//     // psi_m200_terra
//     if (L[col] < 0)
//       psi_200_line[col] = 2 * log((1 + x_200_line[col]) / 2) + log((1 + x_200_line[col] * x_200_line[col]) / 2) - 2 * atan(x_200_line[col]) + 0.5 * PI;
//     else
//       psi_200_line[col] = -5 * ((10 - DISP) / L[col]);

//     // psi_m2_terra
//     if (L[col] < 0)
//       psi_2_line[col] = 2 * log((1 + y_2_line[col] * y_2_line[col]) / 2);
//     else
//       psi_2_line[col] = -5 * ((10 - DISP) / L[col]);

//     // u*
//     double ust = (VON_KARMAN * ustar_previous[line][col]) / (log((10 - DISP) / zom_vector[line][col]) - psi_200_line[col]);

//     // rah
//     double zoh_terra = zom_vector[line][col] / pow(exp(1.0), (kb1_vector[line][col]));
//     double temp_rah1_corr_terra = (ust * VON_KARMAN);
//     double temp_rah2_corr_terra = log((10 - DISP) / zom_vector[line][col]) - psi_2_line[col];
//     double temp_rah3_corr_terra = temp_rah1_corr_terra * log(zom_vector[line][col] / zoh_terra);
//     double rah = (temp_rah1_corr_terra * temp_rah2_corr_terra) + temp_rah3_corr_terra;

//     ustar_vector[line][col] = ust;
//     aerodynamic_resistance_vector[line][col] = rah;

//     if (line == hot_pixel.line && col == hot_pixel.col)
//     {
//       hot_pixel.aerodynamic_resistance.push_back(rah);
//     }

//     if (line == cold_pixel.line && col == cold_pixel.col)
//     {
//       cold_pixel.aerodynamic_resistance.push_back(rah);
//     }
//   }
// }