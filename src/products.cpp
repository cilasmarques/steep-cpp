#include "products.h"

void d0_fuction(std::vector<double> pai_line, int width_band, double d0_line[])
{
  double CD1 = 20.6;
  double HGHT = 4; // TODO: PRECISO DO HGHT (altura da vegetação, acredito que seja um dado de input)!!

  for (int col = 0; col < width_band; col++)
  {

    double pai = pai_line[col];
    double cd1_pai_root = sqrt(CD1 * pai);

    double DISP = HGHT * ((1 - (1 / cd1_pai_root)) + (pow(exp(1.0), -cd1_pai_root) / cd1_pai_root));
    if (pai < 0) {
      DISP = 0;
    }

    d0_line[col] = DISP;
  }
};

void kb_function(double ustar_line[], double zom_line[], std::vector<double> pai_line, double ndvi_line[], double ndvi_max, double ndvi_min, int width_band, double kb1_line[])
{
  double HGHT = 4; // TODO: PRECISO DO HGHT (altura da vegetação, acredito que seja um dado de input)!!

  double visc = 0.00001461;
  double pr = 0.71;
  double c1 = 0.320;
  double c2 = 0.264;
  double c3 = 15.1;
  double cd = 0.2;
  double ct = 0.01;
  double sf_c = 0.3;
  double sf_d = 2.5;
  double sf_e = 4.0;

  for (int col = 0; col < width_band; col++)
  {
    double zom = zom_line[col];
    double u_ast_ini_terra = ustar_line[col];
    double PAI = pai_line[col];

    double Re_star = (u_ast_ini_terra * 0.009) / visc;
    double Ct_star = pow(pr, -0.667) * pow(Re_star, -0.5);
    double beta = c1 - c2 * (exp((cd * -c3 * PAI)));
    double nec_terra = (cd * PAI) / (beta * beta * 2);

    double kb1_fst_part = (cd * VON_KARMAN) / (4 * ct * beta * (1 - exp(nec_terra * -0.5)));
    double kb1_sec_part = (beta * VON_KARMAN * (zom / HGHT)) / Ct_star;
    double kb1s = (pow(Re_star, 0.25) * 2.46) - 2;

    double fc = 1 - pow((ndvi_line[col] - ndvi_max) / (ndvi_min - ndvi_max), 0.4631);
    double fs = 1 - fc;

    double soil_moisture_day_rel = 0.33; // TODO: PRECISO DO soil mpoisture - SM (altura da vegetação, acredito que seja um dado de input)!!

    double SF = sf_c + (1 / (1 + pow(exp(1.0), (sf_d - (sf_e * soil_moisture_day_rel)))));

    kb1_line[col] = ((kb1_fst_part * pow(fc, 2)) +
                     (kb1_sec_part * pow(fc, 2) * pow(fs, 2)) +
                     (pow(fs, 2) * kb1s)) *
                    SF;
  }
};

/**
 * @brief  Computes the momentum roughness length (zom).
 * @param  A_ZOM: Correlation constant a.
 * @param  B_ZOM: Correlation constant b.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  zom_line[]: Auxiliary array for save the calculated value of zom for the line.
 */
void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[], double d0_line[], std::vector<double> pai_line)
{
  double HGHT = 4; // TODO: PRECISO DO HGHT (altura da vegetação, acredito que seja um dado de input)!!
  double CD = 0.01;
  double CR = 0.35;
  double PSICORR = 0.2;
  
  double gama;
  for (int col = 0; col < width_band; col++)
  {
    gama = pow((CD + CR * (pai_line[col] / 2)), -0.5);

    if (gama < 3.3)
      gama = 3.3;

    zom_line[col] = (HGHT - d0_line[col]) * pow(exp(1.0), (-VON_KARMAN * gama) + PSICORR);
  }
};

/**
 * @brief  The friction velocity (u*) is computed.
 * @param  u10: Wind speed at 10 m.
 * @param  zom_line[]: Array containing the specified line from the zom computation.
 * @param  width_band: Band width.
 * @param  ustar_line[]: Auxiliary array for save the calculated value of ustar for the line.
 */
void ustar_fuction(double u10, double zom_line[], int width_band, double ustar_line[], double d0_line[])
{
  double zu = 10; // Altura da velocidade do vento

  for (int col = 0; col < width_band; col++)
  {
    double DISP = d0_line[col]; // Deslocamento rugoso
    double zom = zom_line[col]; // Altura do comprimento de mistura
    ustar_line[col] = (u10 * VON_KARMAN) / std::log((zu - DISP) / zom);
  }
};

/**
 * @brief  Computes the aerodynamic resistance (Rah).
 * @param  ustar_line[]: Array containing the specified line from the ustar computation.
 * @param  width_band: Band width.
 * @param  aerodynamic_resistance_line[]: Auxiliary array for save the calculated value of Rah for the line.
 */
void aerodynamic_resistance_fuction(double ustar_line[], int width_band, double aerodynamic_resistance_line[], double zom_line[], double d0_line[], double kb1_line[])
{
  double zu = 10.0; // Altura da velocidade do vento

  for (int col = 0; col < width_band; col++)
  {
    double DISP = d0_line[col]; // Deslocamento rugoso
    double zom = zom_line[col]; // Altura do comprimento de mistura
    double zoh_terra = zom / pow(exp(1.0), (kb1_line[col]));


    // Codigo STEEP
    double temp_kb_1_terra = log(zom / zoh_terra);
    double temp_rah1_terra = (1 / (ustar_line[col] * VON_KARMAN));
    double temp_rah2 = log(((zu - DISP) / zom));
    double temp_rah3_terra = temp_rah1_terra * temp_kb_1_terra;

    aerodynamic_resistance_line[col] = temp_rah1_terra * temp_rah2 + temp_rah3_terra;

    // Artigo
    // aerodynamic_resistance_line[col] = (1 / (ustar_line[col] * VON_KARMAN)) * log(((zu - DISP) / zom)) + ((1 / (ustar_line[col] * VON_KARMAN)) * kb1_line[col]);
  }
};

/**
 * @brief  Computes Latent Heat Flux (LE).
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_flux_line[]: Array containing the specified line from the G computation.
 * @param  sensible_heat_flux_line[]: Array containing the specified line from the H computation.
 * @param  width_band: Band width.
 * @param  latent_heat_flux[]: Auxiliary array for save the calculated value of LE for the line.
 */
void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[])
{
  for (int col = 0; col < width_band; col++)
    latent_heat_flux[col] = net_radiation_line[col] - soil_heat_flux_line[col] - sensible_heat_flux_line[col];
};

/**
 * @brief  Calculates the Net Radiation for 24 hours (Rn24h).
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  Ra24h: Extraterrestrial Radiation defined as solar short wave radiation in the absence of an atmosphere (Ra24h).
 * @param  Rs24h: Short wave radiation incident in 24 hours (Rs24h).
 * @param  width_band: Band width.
 * @param  net_radiation_24h_line[]: Auxiliary array for save the calculated value of Rn24h for the line.
 */
void net_radiation_24h_function(double albedo_line[], double Ra24h, double Rs24h, int width_band, double net_radiation_24h_line[])
{
  int FL = 110;

  for (int col = 0; col < width_band; col++)
    net_radiation_24h_line[col] = (1 - albedo_line[col]) * Rs24h - FL * Rs24h / Ra24h;
};

/**
 * @brief  The Reference ET Fraction (EF) is computed.
 * @param  latent_heat_flux_line[]: Array containing the specified line from the LE computation.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_line[]: Array containing the specified line from the G computation.
 * @param  width_band: Band width.
 * @param  evapotranspiration_fraction_line[]: Auxiliary array for save the calculated value of EF for the line.
 */
void evapotranspiration_fraction_fuction(double latent_heat_flux_line[], double net_radiation_line[], double soil_heat_line[], int width_band, double evapotranspiration_fraction_line[])
{
  for (int col = 0; col < width_band; col++) {
    double latent = latent_heat_flux_line[col];
    double net_rad = net_radiation_line[col];
    double soil_heat = soil_heat_line[col];

    evapotranspiration_fraction_line[col] = latent / (net_rad - soil_heat);
  }
};

/**
 * @brief  Computes Sensible Heat Flux for 24 hours (H24h).
 * @param  evapotranspiration_fraction_line[]: Array containing the specified line from the EF computation.
 * @param  net_radiation_24h_line[]: Array containing the specified line from the Rn24h computation.
 * @param  width_band: Band width.
 * @param  sensible_heat_flux_24h_line[]: Auxiliary array for save the calculated value of H24h for the line.
 */
void sensible_heat_flux_24h_fuction(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double sensible_heat_flux_24h_line[])
{
  for (int col = 0; col < width_band; col++)
    sensible_heat_flux_24h_line[col] = (1 - evapotranspiration_fraction_line[col]) * net_radiation_24h_line[col];
};

/**
 * @brief  Calculates Latente Heat Flux for 24 hours (LE24h).
 * @param  evapotranspiration_fraction_line[]: Array containing the specified line from the EF computation.
 * @param  net_radiation_24h_line[]: Array containing the specified line from the Rn24h computation.
 * @param  width_band: Band width.
 * @param  latent_heat_flux_24h_line[]: Auxiliary array for save the calculated value of LE24h for the line.
 */
void latent_heat_flux_24h_function(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double latent_heat_flux_24h_line[])
{
  for (int col = 0; col < width_band; col++)
    latent_heat_flux_24h_line[col] = evapotranspiration_fraction_line[col] * net_radiation_24h_line[col];
};

/**
 * @brief  Computes the Evapotranspiration for 24 hours (ET24h)
 * @param  latent_heat_flux_24h_line[]: Array containing the specified line from the LE24h computation.
 * @param  station: Station struct.
 * @param  width_band: Band width.
 * @param  evapotranspiration_24h_line[]: Auxiliary array for save the calculated value of ET24h for the line.
 */
void evapotranspiration_24h_function(double latent_heat_flux_24h_line[], Station station, int width_band, double evapotranspiration_24h_line[])
{
  for (int col = 0; col < width_band; col++)
    evapotranspiration_24h_line[col] = (latent_heat_flux_24h_line[col] * 86400) / ((2.501 - 0.00236 * (station.v7_max + station.v7_min) / 2) * 1e+6);
};


void evapotranspiration_ulisses_function(double net_radiation_24h_line[], double evapotranspiration_fraction_line[], int width_band, double evapotranspiration_24h_line[])
{
  for (int col = 0; col < width_band; col++)
    evapotranspiration_24h_line[col] =  net_radiation_24h_line[col] * evapotranspiration_fraction_line[col] * 0.035;
  
};
