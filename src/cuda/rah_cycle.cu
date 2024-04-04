#include "rah_cycle.cuh"

__global__ void correctionCycle(double *surfaceTemperatureLine, double *zomLine, double *ustarRLine, double *ustarWLine, double *rahRLine,
                                double *rahWLine, double *a, double *b, double *u200, int *size)
{
  // Identify position
  int pos = threadIdx.x + blockIdx.x * blockDim.x;

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

    ustarWLine[pos] = (VON_KARMAN * ustarRLine[pos]) / (log(200 / zomLine[pos]) - psi200);
    rahWLine[pos] = (log(2 / 0.1) - psi2 + psi01) / (ustarWLine[pos] * VON_KARMAN);

    pos += blockDim.x * gridDim.x;
  }
}

__global__ void correctionCycleSTEEP(double *surface_temperature_pointer, double *d0_pointer, double *kb1_pointer, double *zom_pointer, double *ustarR_pointer,
                                     double *ustarW_pointer, double *rahR_pointer, double *rahW_pointer, double *H_pointer, double a, double b, int height,
                                     int width)
{
  // Identify position
  unsigned int col = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int row = threadIdx.y + blockIdx.y * blockDim.y;

  while (row < height)
  {
    if (col < width && row < height)
    {
      unsigned int pos = row * width + col;

      double DISP = d0_pointer[pos];
      double dT_ini_terra = a + b * (surface_temperature_pointer[pos] - 273.15);

      double sensibleHeatFlux = RHO * SPECIFIC_HEAT_AIR * (dT_ini_terra) / rahR_pointer[pos];
      double L = -1 * ((RHO * SPECIFIC_HEAT_AIR * pow(ustarR_pointer[pos], 3) * surface_temperature_pointer[pos]) / (VON_KARMAN * GRAVITY * sensibleHeatFlux));

      double y2 = pow((1 - (16 * (10 - DISP)) / L), 0.25);
      double x200 = pow((1 - (16 * (10 - DISP)) / L), 0.25);

      double psi2, psi200;
      if (!isnan(L) && L > 0)
      {
        psi2 = -5 * ((10 - DISP) / L);
        psi200 = -5 * ((10 - DISP) / L);
      }
      else
      {
        psi2 = 2 * log((1 + y2 * y2) / 2);
        psi200 = 2 * log((1 + x200) / 2) + log((1 + x200 * x200) / 2) - 2 * atan(x200) + 0.5 * M_PI;
      }

      double ust = (VON_KARMAN * ustarR_pointer[pos]) / (log((10 - DISP) / zom_pointer[pos]) - psi200);

      double zoh_terra = zom_pointer[pos] / pow(exp(1.0), (kb1_pointer[pos]));
      double temp_rah1_corr_terra = (ust * VON_KARMAN);
      double temp_rah2_corr_terra = log((10 - DISP) / zom_pointer[pos]) - psi2;
      double temp_rah3_corr_terra = temp_rah1_corr_terra * log(zom_pointer[pos] / zoh_terra);
      double rah = (temp_rah1_corr_terra * temp_rah2_corr_terra) + temp_rah3_corr_terra;

      ustarW_pointer[pos] = ust;
      rahW_pointer[pos] = rah;
      H_pointer[pos] = sensibleHeatFlux;
    }
    row += blockDim.y * gridDim.y;
  }
}