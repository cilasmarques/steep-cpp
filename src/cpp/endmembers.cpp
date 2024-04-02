#include "endmembers.h"

void compute_H0(vector<double> net_radiation_line, vector<double> soil_heat_flux, int width_band, vector<double> &ho_line)
{
  for (int col = 0; col < width_band; col++)
    ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];
};

void calculate_quartiles(double *values, int size, double *quartiles, const vector<double> &percentiles)
{
  for (size_t i = 0; i < percentiles.size(); ++i)
  {
    int nth = int(percentiles[i] * size);
    nth_element(values, values + nth, values + size);
    quartiles[i] = values[nth];
  }
}

void get_quartiles(vector<vector<double>> ndvi_vector, vector<vector<double>> albedo_vector, vector<vector<double>> surface_temperature_vector, int height_band, int width_band, double *ndvi_quartile, double *albedo_quartile, double *ts_quartile)
{
  const int SIZE = height_band * width_band;

  double *ndvi_values = new double[SIZE];
  double *albedo_values = new double[SIZE];
  double *ts_values = new double[SIZE];

  int ndvi_count = 0, albedo_count = 0, ts_count = 0;

  for (int i = 0; i < height_band; ++i)
  {
    for (int j = 0; j < width_band; ++j)
    {
      if (!isnan(ndvi_vector[i][j]) && !isinf(ndvi_vector[i][j]))
      {
        ndvi_values[ndvi_count++] = ndvi_vector[i][j];
      }
      if (!isnan(albedo_vector[i][j]) && !isinf(albedo_vector[i][j]))
      {
        albedo_values[albedo_count++] = albedo_vector[i][j];
      }
      if (!isnan(surface_temperature_vector[i][j]) && !isinf(surface_temperature_vector[i][j]))
      {
        ts_values[ts_count++] = surface_temperature_vector[i][j];
      }
    }
  }

  // Definindo os percentis para calcular os quartis desejados
  vector<double> ndvi_percentiles = {0.15, 0.85, 0.97};
  vector<double> albedo_percentiles = {0.25, 0.50, 0.75};
  vector<double> ts_percentiles = {0.20, 0.85, 0.97};

  calculate_quartiles(ndvi_values, ndvi_count, ndvi_quartile, ndvi_percentiles);
  calculate_quartiles(albedo_values, albedo_count, albedo_quartile, albedo_percentiles);
  calculate_quartiles(ts_values, ts_count, ts_quartile, ts_percentiles);

  delete[] ndvi_values;
  delete[] albedo_values;
  delete[] ts_values;
}

pair<Candidate, Candidate> getColdHotPixelsSTEPP(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band)
{
  vector<Candidate> hotCandidatesGroup;
  vector<Candidate> coldCandidatesGroup;
  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  vector<double> hotNDVIQuartile(3);
  vector<double> hotTSQuartile(3);
  vector<double> hotAlbedoQuartile(3);

  get_quartiles(ndvi_vector, albedo_vector, surface_temperature_vector, height_band, width_band, hotNDVIQuartile.data(), hotAlbedoQuartile.data(), hotTSQuartile.data());

  for (int line = 0; line < height_band; line++)
  {
    vector<double> ho_line = ho_vector[line];
    vector<double> ndvi_line = ndvi_vector[line];
    vector<double> surface_temperature_line = surface_temperature_vector[line];
    vector<double> albedo_line = albedo_vector[line];
    vector<double> net_radiation_line = net_radiation_vector[line];
    vector<double> soil_heat_line = soil_heat_vector[line];

    compute_H0(net_radiation_line, soil_heat_line, width_band, ho_line);

    for (int col = 0; col < width_band; col++)
    {
      bool ndviValid = !std::isnan(ndvi_line[col]) && ndvi_line[col] > 0.10 && ndvi_line[col] < hotNDVIQuartile[0];
      bool albedoValid = !std::isnan(albedo_line[col]) && albedo_line[col] > hotAlbedoQuartile[1] && albedo_line[col] < hotAlbedoQuartile[2];
      bool tsValid = !std::isnan(surface_temperature_line[col]) && surface_temperature_line[col] > hotTSQuartile[1] && surface_temperature_line[col] < hotTSQuartile[2];

      if (albedoValid && ndviValid && tsValid)
        hotCandidatesGroup.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);

      ndviValid = !std::isnan(ndvi_line[col]) && ndvi_line[col] > hotNDVIQuartile[2];
      albedoValid = !std::isnan(surface_temperature_line[col]) && albedo_line[col] > hotAlbedoQuartile[0] && albedo_line[col] < hotAlbedoQuartile[1];
      tsValid = !std::isnan(albedo_line[col]) && surface_temperature_line[col] < hotTSQuartile[0];

      if (albedoValid && ndviValid && tsValid)
        coldCandidatesGroup.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (hotCandidatesGroup.empty() || coldCandidatesGroup.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  std::sort(hotCandidatesGroup.begin(), hotCandidatesGroup.end(), compare_candidate_temperature);
  unsigned int posHot = static_cast<unsigned int>(std::floor(hotCandidatesGroup.size() * 0.5));

  std::sort(coldCandidatesGroup.begin(), coldCandidatesGroup.end(), compare_candidate_temperature);
  unsigned int posCold = static_cast<unsigned int>(std::floor(coldCandidatesGroup.size() * 0.5));

  return {hotCandidatesGroup[posHot], coldCandidatesGroup[posCold]};
}
