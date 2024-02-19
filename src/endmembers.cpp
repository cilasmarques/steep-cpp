#include "endmembers.h"

void compute_H0(vector<double> net_radiation_line, vector<double> soil_heat_flux, int width_band, vector<double> &ho_line)
{
  for (int col = 0; col < width_band; col++)
    ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];
};

void filter_valid_values(vector<double> target_line, double *target_values, int width_band, int *pos)
{
  for (int col = 0; col < width_band; col++)
  {
    if (!isnan(target_line[col]) && !isinf(target_line[col]))
    {
      target_values[*pos] = target_line[col];
      (*pos)++;
    }
  }
}

void get_quartiles(vector<vector<double>> target_vector, double *v_quartile, int height_band, int width_band, double first_interval, double mid_interval, double last_interval)
{
  const int SIZE = height_band * width_band;
  double *target_values = (double *)malloc(sizeof(double) * SIZE);

  if (target_values == NULL)
    exit(15);

  int pos = 0;
  for (int line = 0; line < height_band; line++)
  {
    filter_valid_values(target_vector[line], target_values, width_band, &pos);
  }

  sort(target_values, target_values + pos);

  v_quartile[0] = target_values[int(floor(first_interval * pos))];
  v_quartile[1] = target_values[int(floor(mid_interval * pos))];
  v_quartile[2] = target_values[int(floor(last_interval * pos))];

  free(target_values);
}

pair<Candidate, Candidate> getColdHotPixelsSTEPP(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band)
{
  vector<Candidate> hotCandidatesGroup;
  vector<Candidate> coldCandidatesGroup;

  vector<double> hotNDVIQuartile(3);
  vector<double> hotTSQuartile(3);
  vector<double> hotAlbedoQuartile(3);

  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  get_quartiles(ndvi_vector, hotNDVIQuartile.data(), height_band, width_band, 0.15, 0.85, 0.97);
  get_quartiles(albedo_vector, hotAlbedoQuartile.data(), height_band, width_band, 0.25, 0.50, 0.75);
  get_quartiles(surface_temperature_vector, hotTSQuartile.data(), height_band, width_band, 0.20, 0.85, 0.97);

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
