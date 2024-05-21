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

void get_quartiles(vector<vector<double>> target_vector, double *v_quartile, int height_band, int width_band, double first_interval, double middle_interval, double last_interval)
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

  int first_index = static_cast<int>(floor(first_interval * pos));
  int middle_index = static_cast<int>(floor(middle_interval * pos));
  int last_index = static_cast<int>(floor(last_interval * pos));

  std::nth_element(target_values, target_values + first_index, target_values + pos);
  v_quartile[0] = target_values[first_index];

  std::nth_element(target_values, target_values + middle_index, target_values + pos);
  v_quartile[1] = target_values[middle_index];

  std::nth_element(target_values, target_values + last_index, target_values + pos);
  v_quartile[2] = target_values[last_index];

  free(target_values);
}


pair<Candidate, Candidate> getEndmembersASEBAL(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band)
{
  vector<Candidate> hotCandidates;
  vector<Candidate> coldCandidates;
  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  vector<double> tsQuartile(3);
  vector<double> ndviQuartile(3);
  vector<double> albedoQuartile(3);
  get_quartiles(ndvi_vector, ndviQuartile.data(), height_band, width_band, 0.25, 0.75, 0.75);
  get_quartiles(albedo_vector, albedoQuartile.data(), height_band, width_band, 0.25, 0.50, 0.75);
  get_quartiles(surface_temperature_vector, tsQuartile.data(), height_band, width_band, 0.25, 0.75, 0.75);

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

      bool hotAlbedo = !isnan(albedo_line[col]) && albedo_line[col] > albedoQuartile[1];
      bool hotNDVI = !isnan(ndvi_line[col]) && ndvi_line[col] > 0.10 && ndvi_line[col] < ndviQuartile[0];
      bool hotTS = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] > tsQuartile[1];

      bool coldAlbedo = !isnan(albedo_line[col]) && albedo_line[col] < albedoQuartile[1];
      bool coldNDVI = !isnan(ndvi_line[col]) && ndvi_line[col] >= ndviQuartile[1];
      bool coldTS = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] < tsQuartile[0];

      if (hotAlbedo && hotNDVI && hotTS)
        hotCandidates.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);

      if (coldAlbedo && coldNDVI && coldTS)
        coldCandidates.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (hotCandidates.empty() || coldCandidates.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  // Creating second pixel group, all values lower than the 3rd quartile are excluded
  std::sort(hotCandidates.begin(), hotCandidates.end(), compare_candidate_temperature);
  std::sort(coldCandidates.begin(), coldCandidates.end(), compare_candidate_temperature);
  
  unsigned int hotPos = static_cast<unsigned int>(std::floor(hotCandidates.size() * 0.5));
  unsigned int coldPos = static_cast<unsigned int>(std::floor(coldCandidates.size() * 0.5));

  return {hotCandidates[hotPos], coldCandidates[coldPos]};
}

pair<Candidate, Candidate> getEndmembersSTEPP(vector<vector<double>> ndvi_vector, vector<vector<double>> surface_temperature_vector, vector<vector<double>> albedo_vector, vector<vector<double>> net_radiation_vector, vector<vector<double>> soil_heat_vector, int height_band, int width_band)
{
  vector<Candidate> hotCandidates;
  vector<Candidate> coldCandidates;

  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  vector<double> tsQuartile(3);
  vector<double> ndviQuartile(3);
  vector<double> albedoQuartile(3);
  get_quartiles(ndvi_vector, ndviQuartile.data(), height_band, width_band, 0.15, 0.85, 0.97);
  get_quartiles(albedo_vector, albedoQuartile.data(), height_band, width_band, 0.25, 0.50, 0.75);
  get_quartiles(surface_temperature_vector, tsQuartile.data(), height_band, width_band, 0.20, 0.85, 0.97);

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
      bool hotNDVI = !std::isnan(ndvi_line[col]) && ndvi_line[col] > 0.10 && ndvi_line[col] < ndviQuartile[0];
      bool hotAlbedo = !std::isnan(albedo_line[col]) && albedo_line[col] > albedoQuartile[1] && albedo_line[col] < albedoQuartile[2];
      bool hotTS = !std::isnan(surface_temperature_line[col]) && surface_temperature_line[col] > tsQuartile[1] && surface_temperature_line[col] < tsQuartile[2];

      bool coldNDVI = !std::isnan(ndvi_line[col]) && ndvi_line[col] > ndviQuartile[2];
      bool coldAlbedo = !std::isnan(surface_temperature_line[col]) && albedo_line[col] > albedoQuartile[0] && albedo_line[col] < albedoQuartile[1];
      bool coldTS = !std::isnan(albedo_line[col]) && surface_temperature_line[col] < tsQuartile[0];

      if (hotAlbedo && hotNDVI && hotTS)
        hotCandidates.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);

      if (coldNDVI && coldAlbedo && coldTS)
        coldCandidates.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (hotCandidates.empty() || coldCandidates.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  std::sort(hotCandidates.begin(), hotCandidates.end(), compare_candidate_temperature);
  std::sort(coldCandidates.begin(), coldCandidates.end(), compare_candidate_temperature);

  unsigned int hotPos = static_cast<unsigned int>(std::floor(hotCandidates.size() * 0.5));
  unsigned int coldPos = static_cast<unsigned int>(std::floor(coldCandidates.size() * 0.5));

  return {hotCandidates[hotPos], coldCandidates[coldPos]};
}
