#include "kernels.h"

void rah_correction_cycle_STEEP(int start_line, int end_line, int width_band, double a, double b, double *surface_temperature_pointer,
                                double *d0_pointer, double *zom_pointer, double *kb1_pointer, double *sensible_heat_flux_pointer, 
                                double *ustar_pointer, double *aerodynamic_resistance_pointer)
{
  for (int line = start_line; line < end_line; line++)
  {
    for (int col = 0; col < width_band; col++)
    {
      int pos = line * width_band + col;

      double DISP = d0_pointer[pos];
      double dT_ini_terra = (a + b * (surface_temperature_pointer[pos] - 273.15));

      // H_ini_terra
      sensible_heat_flux_pointer[pos] = RHO * SPECIFIC_HEAT_AIR * (dT_ini_terra) / aerodynamic_resistance_pointer[pos];

      // L_MB_terra
      double ustar_pow_3 = ustar_pointer[pos] * ustar_pointer[pos] * ustar_pointer[pos];

      double L = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_pointer[pos]) / (VON_KARMAN * GRAVITY * sensible_heat_flux_pointer[pos]));

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

      // u*
      double ust = (VON_KARMAN * ustar_pointer[pos]) / (log((10 - DISP) / zom_pointer[pos]) - psi200);

      // rah
      double zoh_terra = zom_pointer[pos] / pow(exp(1.0), (kb1_pointer[pos]));
      double temp_rah1_corr_terra = (ust * VON_KARMAN);
      double temp_rah2_corr_terra = log((10 - DISP) / zom_pointer[pos]) - psi2;
      double temp_rah3_corr_terra = temp_rah1_corr_terra * log(zom_pointer[pos] / zoh_terra);
      double rah = (temp_rah1_corr_terra * temp_rah2_corr_terra) + temp_rah3_corr_terra;

      ustar_pointer[pos] = ust;
      aerodynamic_resistance_pointer[pos] = rah;
    }
  }
};
