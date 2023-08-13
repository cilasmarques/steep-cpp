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

void get_quartiles(vector<vector<double>> target_matriz, double *v_quartile, int height_band, int width_band, double first_interval, double last_interval)
{
  const int SIZE = height_band * width_band;
  double *target_values = (double *)malloc(sizeof(double) * SIZE);

  if (target_values == NULL)
    exit(15);

  int pos = 0;
  for (int line = 0; line < height_band; line++)
  {
    filter_valid_values(target_matriz[line], target_values, width_band, &pos);
  }

  sort(target_values, target_values + pos);

  v_quartile[0] = target_values[int(floor(first_interval * pos))];
  v_quartile[1] = target_values[int(floor(last_interval * pos))];

  free(target_values);
}

Candidate getHotPixelASEBAL(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band)
{
  vector<Candidate> candidatesGroup;
  vector<double> ndviQuartile(2);
  vector<double> tsQuartile(2);
  vector<double> albedoQuartile(2);
  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  get_quartiles(ndvi_matriz, ndviQuartile.data(), height_band, width_band, 0.25, 0.75);
  get_quartiles(albedo_matriz, albedoQuartile.data(), height_band, width_band, 0.25, 0.75);
  get_quartiles(surface_temperature_matriz, tsQuartile.data(), height_band, width_band, 0.25, 0.75);

  for (int line = 0; line < height_band; line++)
  {
    vector<double> ho_line = ho_vector[line];
    vector<double> ndvi_line = ndvi_matriz[line];
    vector<double> surface_temperature_line = surface_temperature_matriz[line];
    vector<double> albedo_line = albedo_matriz[line];
    vector<double> net_radiation_line = net_radiation_matriz[line];
    vector<double> soil_heat_line = soil_heat_matriz[line];

    compute_H0(net_radiation_line, soil_heat_line, width_band, ho_line);

    for (int col = 0; col < width_band; col++)
    {

      bool albedoValid = !isnan(albedo_line[col]) && albedo_line[col] > albedoQuartile[1];
      bool ndviValid = !isnan(ndvi_line[col]) && ndvi_line[col] > 0.10 && ndvi_line[col] < ndviQuartile[0];
      bool tsValid = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] > tsQuartile[1];

      if (albedoValid && ndviValid && tsValid)
        candidatesGroup.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (candidatesGroup.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  // Creating second pixel group, all values lower than the 3rd quartile are excluded
  std::sort(candidatesGroup.begin(), candidatesGroup.end(), compare_candidate_temperature);
  unsigned int pos = static_cast<unsigned int>(std::floor(candidatesGroup.size() * 0.75));
  vector<Candidate> candidatesFinalGroup(candidatesGroup.begin() + pos, candidatesGroup.end());

  if (candidatesFinalGroup.size() <= 0)
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  pos = int(floor(candidatesFinalGroup.size() * 0.5));
  Candidate hotPixel = candidatesFinalGroup[pos];

  return hotPixel;
}

Candidate getColdPixelASEBAL(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band)
{
  vector<Candidate> candidatesGroup;
  vector<double> ndviQuartile(2);
  vector<double> tsQuartile(2);
  vector<double> albedoQuartile(2);
  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  get_quartiles(ndvi_matriz, ndviQuartile.data(), height_band, width_band, 0.25, 0.75);
  get_quartiles(albedo_matriz, albedoQuartile.data(), height_band, width_band, 0.25, 0.50);
  get_quartiles(surface_temperature_matriz, tsQuartile.data(), height_band, width_band, 0.25, 0.75);

  for (int line = 0; line < height_band; line++)
  {
    vector<double> ho_line = ho_vector[line];
    vector<double> ndvi_line = ndvi_matriz[line];
    vector<double> surface_temperature_line = surface_temperature_matriz[line];
    vector<double> albedo_line = albedo_matriz[line];
    vector<double> net_radiation_line = net_radiation_matriz[line];
    vector<double> soil_heat_line = soil_heat_matriz[line];

    compute_H0(net_radiation_line, soil_heat_line, width_band, ho_line);

    for (int col = 0; col < width_band; col++)
    {

      bool albedoValid = !isnan(albedo_line[col]) && albedo_line[col] < albedoQuartile[1];
      bool ndviValid = !isnan(ndvi_line[col]) && ndvi_line[col] >= ndviQuartile[1]; // ndvi_line[col] >= ndviQuartile[3];
      bool tsValid = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] < tsQuartile[0];

      if (albedoValid && ndviValid && tsValid)
        candidatesGroup.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (candidatesGroup.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  // Creating second pixel group, all values lower than the 3rd quartile are excluded
  std::sort(candidatesGroup.begin(), candidatesGroup.end(), compare_candidate_temperature);
  unsigned int pos = static_cast<unsigned int>(std::floor(candidatesGroup.size() * 0.25));
  vector<Candidate> candidatesFinalGroup(candidatesGroup.begin(), candidatesGroup.end() + pos);

  if (candidatesFinalGroup.size() <= 0)
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  pos = int(floor(candidatesFinalGroup.size() * 0.5));
  Candidate coldPixel = candidatesFinalGroup[pos];

  return coldPixel;
}

Candidate getHotPixelSTEPP(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band)
{
  vector<Candidate> candidatesGroup;
  vector<double> ndviQuartile(2);
  vector<double> tsQuartile(2);
  vector<double> albedoQuartile(2);
  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  get_quartiles(ndvi_matriz, ndviQuartile.data(), height_band, width_band, 0.15, 0.85);
  get_quartiles(albedo_matriz, albedoQuartile.data(), height_band, width_band, 0.50, 0.75);
  get_quartiles(surface_temperature_matriz, tsQuartile.data(), height_band, width_band, 0.85, 0.97);

  for (int line = 0; line < height_band; line++)
  {
    vector<double> ho_line = ho_vector[line];
    vector<double> ndvi_line = ndvi_matriz[line];
    vector<double> surface_temperature_line = surface_temperature_matriz[line];
    vector<double> albedo_line = albedo_matriz[line];
    vector<double> net_radiation_line = net_radiation_matriz[line];
    vector<double> soil_heat_line = soil_heat_matriz[line];

    compute_H0(net_radiation_line, soil_heat_line, width_band, ho_line);

    for (int col = 0; col < width_band; col++)
    {
      bool ndviValid = !std::isnan(ndvi_line[col]) && ndvi_line[col] > 0.10 && ndvi_line[col] < ndviQuartile[0];
      bool albedoValid = !std::isnan(albedo_line[col]) && albedo_line[col] > albedoQuartile[0] && albedo_line[col] < albedoQuartile[1];
      bool tsValid = !std::isnan(surface_temperature_line[col]) && surface_temperature_line[col] > tsQuartile[0] && surface_temperature_line[col] < tsQuartile[1];

      if (albedoValid && ndviValid && tsValid)
        candidatesGroup.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (candidatesGroup.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  std::sort(candidatesGroup.begin(), candidatesGroup.end(), compare_candidate_temperature);
  unsigned int pos = static_cast<unsigned int>(std::floor(candidatesGroup.size() * 0.5));

  return candidatesGroup[pos];
}

Candidate getColdPixelSTEPP(vector<vector<double>> ndvi_matriz, vector<vector<double>> surface_temperature_matriz, vector<vector<double>> albedo_matriz, vector<vector<double>> net_radiation_matriz, vector<vector<double>> soil_heat_matriz, int height_band, int width_band)
{
  vector<Candidate> candidatesGroup;
  vector<double> ndviQuartile(2);
  vector<double> tsQuartile(2);
  vector<double> albedoQuartile(2);
  vector<vector<double>> ho_vector(height_band, vector<double>(width_band));

  get_quartiles(ndvi_matriz, ndviQuartile.data(), height_band, width_band, 0.15, 0.97);
  get_quartiles(albedo_matriz, albedoQuartile.data(), height_band, width_band, 0.25, 0.50);
  get_quartiles(surface_temperature_matriz, tsQuartile.data(), height_band, width_band, 0.20, 0.85);

  for (int line = 0; line < height_band; line++)
  {
    vector<double> ho_line = ho_vector[line];
    vector<double> ndvi_line = ndvi_matriz[line];
    vector<double> surface_temperature_line = surface_temperature_matriz[line];
    vector<double> albedo_line = albedo_matriz[line];
    vector<double> net_radiation_line = net_radiation_matriz[line];
    vector<double> soil_heat_line = soil_heat_matriz[line];

    compute_H0(net_radiation_line, soil_heat_line, width_band, ho_line);

    for (int col = 0; col < width_band; col++)
    {

      bool ndviValid = !std::isnan(ndvi_line[col]) && ndvi_line[col] > ndviQuartile[1];
      bool albedoValid = !std::isnan(surface_temperature_line[col]) && albedo_line[col] > albedoQuartile[0] && albedo_line[col] < albedoQuartile[1];
      bool tsValid = !std::isnan(albedo_line[col]) && surface_temperature_line[col] < tsQuartile[0];

      if (albedoValid && ndviValid && tsValid)
        candidatesGroup.emplace_back(ndvi_line[col], surface_temperature_line[col], net_radiation_line[col], soil_heat_line[col], ho_line[col], line, col);
    }
  }

  if (candidatesGroup.empty())
  {
    cerr << "Pixel problem! - There are no final candidates";
    exit(15);
  }

  std::sort(candidatesGroup.begin(), candidatesGroup.end(), compare_candidate_temperature);
  unsigned int pos = static_cast<unsigned int>(std::floor(candidatesGroup.size() * 0.5));

  return candidatesGroup[pos];
}

// pair<Candidate, Candidate> getColdHotPixelsESA(TIFF **ndvi, TIFF **surface_temperature, TIFF **albedo, TIFF **net_radiation, TIFF **soil_heat, TIFF **landCover, int height_band, int width_band, string output_path)
// {

//   // Testing land cover homogeneity
//   TIFF *outputLC = TIFFOpen((output_path + "/outLC.tif").c_str(), "w8m");
//   setup(outputLC, width_band, height_band, 32, 2);
//   testLandCoverHomogeneity(*landCover, outputLC);
//   TIFFClose(outputLC);

//   // Testing the other params homogeneity
//   outputLC = TIFFOpen((output_path + "/outLC.tif").c_str(), "r");
//   TIFF *outputH = TIFFOpen((output_path + "/outH.tif").c_str(), "w8m");
//   setup(outputH, width_band, height_band, 32, 2);
//   testHomogeneity(*ndvi, *surface_temperature, *albedo, outputLC, outputH);
//   TIFFClose(outputLC);
//   TIFFClose(outputH);

//   // Morphological test
//   outputH = TIFFOpen((output_path + "/outH.tif").c_str(), "r");
//   TIFF *outputAll = TIFFOpen((output_path + "/outAll.tif").c_str(), "w8m");
//   setup(outputAll, width_band, height_band, 32, 2);
//   testMorphological(outputH, outputAll, 50);
//   TIFFClose(outputH);
//   TIFFClose(outputAll);

//   // Getting the candidates
//   outputAll = TIFFOpen((output_path + "/outAll.tif").c_str(), "r");
//   int all_condition[width_band];
//   vector<Candidate> listTS;

//   // Auxiliary arrays
//   double ndvi_line[width_band], surface_temperature_line[width_band];
//   double net_radiation_line[width_band], soil_heat_line[width_band];
//   double ho_line[width_band];

//   // Creating candidates array for TS and then for NDVI as a copy
//   for (int line = 0; line < height_band; line++)
//   {

//     read_line_tiff(outputAll, all_condition, line);

//     read_line_tiff(*net_radiation, net_radiation_line, line);
//     read_line_tiff(*soil_heat, soil_heat_line, line);

//     hoCalc(net_radiation_line, soil_heat_line, width_band, ho_line);

//     read_line_tiff(*ndvi, ndvi_line, line);
//     read_line_tiff(*surface_temperature, surface_temperature_line, line);

//     for (int col = 0; col < width_band; col++)
//     {

//       if (all_condition[col] && !isnan(ndvi_line[col]))
//       {

//         listTS.push_back(Candidate(ndvi_line[col],
//                                    surface_temperature_line[col],
//                                    net_radiation_line[col],
//                                    soil_heat_line[col],
//                                    ho_line[col],
//                                    line, col));
//       }
//     }
//   }

//   if (listTS.size() <= 0)
//   {
//     cerr << "Pixel problem! - There are no precandidates";
//     exit(15);
//   }

//   vector<Candidate> listNDVI(listTS);

//   sort(listTS.begin(), listTS.end(), compare_candidate_temperature);
//   sort(listNDVI.begin(), listNDVI.end(), compare_candidate_ndvi);

//   double ts_min = listTS[0].temperature, ts_max = listTS[listTS.size() - 1].temperature;
//   double ndvi_min = listNDVI[0].ndvi, ndvi_max = listNDVI[listNDVI.size() - 1].ndvi;
//   int binTS = int(ceil((ts_max - ts_min) / 0.25));       // 0.25 is TS bin size
//   int binNDVI = int(ceil((ndvi_max - ndvi_min) / 0.01)); // 0.01 is ndvi bin size

//   vector<Candidate> histTS[binTS], final_histTS;

//   for (Candidate c : listTS)
//   {

//     int pos = int(ceil((c.temperature - ts_min) / 0.25));
//     histTS[pos > 0 ? pos - 1 : 0].push_back(c);
//   }

//   for (int i = 0; i < binTS; i++)
//   {

//     if (histTS[i].size() > 50)
//     {

//       for (Candidate c : histTS[i])
//         final_histTS.push_back(c);
//     }
//   }

//   if (final_histTS.size() <= 0)
//   {
//     cerr << "Pixel problem! - There are no final TS candidates";
//     exit(15);
//   }

//   vector<Candidate> histNDVI[binNDVI], final_histNDVI;
//   for (Candidate c : listNDVI)
//   {

//     int pos = int(ceil((c.ndvi - ndvi_min) / 0.01));
//     histNDVI[pos > 0 ? pos - 1 : 0].push_back(c);
//   }

//   for (int i = 0; i < binNDVI; i++)
//   {

//     if (histNDVI[i].size() > 50)
//     {

//       for (Candidate c : histNDVI[i])
//         final_histNDVI.push_back(c);
//     }
//   }

//   if (final_histNDVI.size() <= 0)
//   {
//     cerr << "Pixel problem! - There are no final NDVI candidates";
//     exit(15);
//   }

//   // Select cold pixel
//   int pixel_count = 0, n1 = 1, n2 = 1, ts_pos, ndvi_pos, beginTs = 0, beginNDVI = final_histNDVI.size() - 1;
//   vector<Candidate> coldPixels;
//   while (pixel_count < 10 && !(n2 == 10 && n1 == 10))
//   {

//     ts_pos = int(floor(n1 / 100.0 * final_histTS.size()));
//     ndvi_pos = int(floor((100 - n2) / 100.0 * final_histNDVI.size()));

//     for (int i = beginTs; i <= ts_pos && pixel_count < 10; i++)
//     {

//       for (int j = beginNDVI; j >= ndvi_pos && pixel_count < 10; j--)
//       {

//         if (equals(final_histTS[i], final_histNDVI[j]))
//         {

//           coldPixels.push_back(final_histTS[i]);
//           pixel_count++;
//         }
//       }
//     }

//     beginTs = ts_pos;
//     beginNDVI = ndvi_pos;

//     if (n2 < 10)
//       n2++;
//     else if (n1 < 10)
//     {
//       n1++;
//       beginNDVI = final_histNDVI.size() - 1;
//     }
//   }

//   if (coldPixels.size() <= 0)
//   {
//     cerr << "Pixel problem! - There are no cold candidates";
//     exit(15);
//   }

//   // Select hot pixel
//   pixel_count = 0, n1 = 1, n2 = 1;
//   vector<Candidate> hotPixels;
//   beginTs = final_histTS.size() - 1, beginNDVI = 0;
//   while (pixel_count < 10 && !(n2 == 10 && n1 == 10))
//   {

//     ts_pos = int(floor((100 - n1) / 100.0 * final_histTS.size()));
//     ndvi_pos = int(floor(n2 / 100.0 * final_histNDVI.size()));

//     for (int i = beginNDVI; i <= ndvi_pos && pixel_count < 10; i++)
//     {

//       for (int j = beginTs; j >= ts_pos && pixel_count < 10; j--)
//       {

//         if (equals(final_histTS[j], final_histNDVI[i]))
//         {

//           hotPixels.push_back(final_histTS[j]);
//           pixel_count++;
//         }
//       }
//     }

//     beginTs = ts_pos;
//     beginNDVI = ndvi_pos;

//     if (n2 < 10)
//       n2++;
//     else if (n1 < 10)
//     {
//       n1++;
//       beginTs = final_histTS.size() - 1;
//     }
//   }

//   if (hotPixels.size() <= 0)
//   {
//     cerr << "Pixel problem! - There are no hot candidates";
//     exit(15);
//   }

//   sort(coldPixels.begin(), coldPixels.end(), compare_candidate_ndvi);
//   sort(hotPixels.begin(), hotPixels.end(), compare_candidate_temperature);

//   return {hotPixels[hotPixels.size() - 1], coldPixels[coldPixels.size() - 1]};
// }
