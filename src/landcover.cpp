// #include "landcover.h"

// /**
//  * @brief   Checks if the input value matches one of the agricultural land cover code.
//  * @param   value: The input value.
//  * @retval  TRUE if the values matches and FALSE otherwise.
//  */
// bool checkLandCode(int value)
// {

//   return (value == AGP) || (value == PAS) || (value == AGR) || (value == CAP) || (value == CSP) || (value == MAP);
// }

// /**
//  * @brief   Tests the land cover homogeneity of a pixel.
//  * @note    A the land cover of a pixel is homogeneous if every neighbour pixel inside a 7x7 window is also agricultural field.
//  * @param   landCover: Land Cover TIFF.
//  * @param   mask: Output binary TIFF, where pixels with 1 means that is a homogeneous pixel and 0 means the otherwise.
//  */
// void testLandCoverHomogeneity(TIFF *landCover, TIFF *mask)
// {

//   uint32 height_band, width_band;
//   TIFFGetField(landCover, TIFFTAG_IMAGELENGTH, &height_band);
//   TIFFGetField(landCover, TIFFTAG_IMAGEWIDTH, &width_band);

//   double **buffer = (double **)malloc(7 * sizeof(double *));

//   for (int i = 0; i < 7; i++)
//   {
//     buffer[i] = (double *)malloc(width_band * sizeof(double));
//   }

//   int relation[7] = {-1, -1, -1, -1, -1, -1, -1}, aux;

//   for (int line = 0; line < height_band; line++)
//   {

//     // Create the respective line of the binary map of eligibles pixels
//     int mask_line[width_band];

//     for (int column = 0; column < width_band; column++)
//     {

//       int pixel_value;

//       aux = line % 7;

//       if (relation[aux] != line)
//       {

//         read_line_tiff(landCover, buffer[aux], line);
//         relation[aux] = line;
//       }

//       pixel_value = buffer[aux][column];

//       mask_line[column] = false;

//       if (checkLandCode(pixel_value))
//       { // Verify if the pixel is an AGR pixel

//         mask_line[column] = true;

//         for (int i = -3; i <= 3 && mask_line[column]; i++)
//         {

//           for (int j = -3; j <= 3 && mask_line[column]; j++)
//           {

//             // Check if the neighbor is AGR too

//             if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band)
//             {

//               aux = (line + j) % 7;

//               if (relation[aux] != (line + j))
//               {

//                 read_line_tiff(landCover, buffer[aux], line + j);
//                 relation[aux] = (line + j);
//               }

//               pixel_value = buffer[aux][column + i];

//               if (!isnan(pixel_value))
//                 if (!checkLandCode(pixel_value))
//                   mask_line[column] = false;
//             }
//           }
//         }
//       }
//     }

//     write_line_tiff(mask, mask_line, line);
//   }

//   for (int i = 0; i < 7; i++)
//   {
//     free(buffer[i]);
//   }
//   free(buffer);
// }

// /**
//  * @brief   Tests the ndvi, surface_temperature and albedo homogeneity of a pixel.
//  * @note    A pixel is homogeneous in these criteria if inside a 7x7 window the coefficient of variation of the albedo and ndvi is less or equal than 25%
//  *          and the surface temperature has a standard deviation less or equal than 1.5 K.
//  * @param   ndvi: NDVI TIFF.
//  * @param   surface_temperature: TS TIFF.
//  * @param   albedo: Albedo TIFF.
//  * @param   maskLC: A binary TIFF conteining the data of the land cover homogeneity.
//  * @param   output: A binary TIFF, where pixels with 1 means that is a homogeneous pixel in land cover, ndvi, surface temperature and albedo, and 0 means otherwise.
//  */
// void testHomogeneity(TIFF *ndvi, TIFF *surface_temperature, TIFF *albedo, TIFF *maskLC, TIFF *output)
// {

//   uint32 height_band, width_band;
//   TIFFGetField(ndvi, TIFFTAG_IMAGELENGTH, &height_band);
//   TIFFGetField(ndvi, TIFFTAG_IMAGEWIDTH, &width_band);

//   double **bufferTS = (double **)malloc(10 * sizeof(double *));
//   double **bufferNDVI = (double **)malloc(10 * sizeof(double *));
//   double **bufferAlb = (double **)malloc(10 * sizeof(double *));

//   for (int i = 0; i < 7; i++)
//   {
//     bufferTS[i] = (double *)malloc(width_band * sizeof(double));
//     bufferNDVI[i] = (double *)malloc(width_band * sizeof(double));
//     bufferAlb[i] = (double *)malloc(width_band * sizeof(double));
//   }

//   int relation[7] = {-1, -1, -1, -1, -1, -1, -1}, aux;

//   for (int line = 0; line < height_band; line++)
//   {

//     // Create the respective line of the binary map of eligibles pixels
//     int mask_line[width_band];
//     read_line_tiff(maskLC, mask_line, line);

//     for (int column = 0; column < width_band; column++)
//     {

//       if (mask_line[column] == true)
//       { // Verify if the pixel passed the land cover test

//         vector<double> ndvi_neighbors;
//         vector<double> ts_neighbors;
//         vector<double> albedo_neighbors;

//         double pixel_value;

//         aux = line % 7;

//         if (relation[aux] != line)
//         {

//           read_line_tiff(surface_temperature, bufferTS[aux], line);
//           read_line_tiff(ndvi, bufferNDVI[aux], line);
//           read_line_tiff(albedo, bufferAlb[aux], line);
//           relation[aux] = line;
//         }

//         if (!isnan(bufferNDVI[aux][column]))
//         {

//           for (int i = -3; i <= 3; i++)
//           {

//             for (int j = -3; j <= 3; j++)
//             {

//               // Add for the NDVI, TS and Albedo the value of neighbors pixels into the respective vector

//               if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band)
//               {

//                 aux = (line + j) % 7;

//                 if (relation[aux] != (line + j))
//                 {

//                   read_line_tiff(surface_temperature, bufferTS[aux], line + j);
//                   read_line_tiff(ndvi, bufferNDVI[aux], line + j);
//                   read_line_tiff(albedo, bufferAlb[aux], line + j);
//                   relation[aux] = line + j;
//                 }

//                 pixel_value = bufferNDVI[aux][column + i];
//                 if (!isnan(pixel_value))
//                   ndvi_neighbors.push_back(pixel_value);

//                 pixel_value = bufferTS[aux][column + i];
//                 if (!isnan(pixel_value))
//                   ts_neighbors.push_back(pixel_value);

//                 pixel_value = bufferAlb[aux][column + i];
//                 if (!isnan(pixel_value))
//                   albedo_neighbors.push_back(pixel_value);
//               }
//             }
//           }

//           // Do the calculation of the dispersion measures from the NDVI, TS and Albedo

//           double meanNDVI, meanTS, meanAlb;
//           double sdNDVI, sdTS, sdAlb;
//           double cvNDVI, cvAlb;
//           double sumNDVI = 0, sumTS = 0, sumAlb = 0;

//           for (int i = 0; i < ndvi_neighbors.size(); i++)
//           {

//             sumNDVI += ndvi_neighbors[i];
//             sumTS += ts_neighbors[i];
//             sumAlb += albedo_neighbors[i];
//           }

//           meanNDVI = sumNDVI / ndvi_neighbors.size();
//           meanTS = sumTS / ts_neighbors.size();
//           meanAlb = sumAlb / albedo_neighbors.size();

//           sumNDVI = 0, sumTS = 0, sumAlb = 0;

//           for (int i = 0; i < ndvi_neighbors.size(); i++)
//           {

//             sumNDVI += (ndvi_neighbors[i] - meanNDVI) * (ndvi_neighbors[i] - meanNDVI);
//             sumTS += (ts_neighbors[i] - meanTS) * (ts_neighbors[i] - meanTS);
//             sumAlb += (albedo_neighbors[i] - meanAlb) * (albedo_neighbors[i] - meanAlb);
//           }

//           sdNDVI = sqrt(sumNDVI / ndvi_neighbors.size());
//           sdTS = sqrt(sumTS / ts_neighbors.size());
//           sdAlb = sqrt(sumAlb / albedo_neighbors.size());

//           cvNDVI = sdNDVI / meanNDVI;
//           cvAlb = sdAlb / meanAlb;

//           // Check if the pixel is eligible
//           mask_line[column] = (cvNDVI < 0.25) && (cvAlb < 0.25) && (sdTS < 1.5);
//         }
//         else
//         {

//           mask_line[column] = false;
//         }
//       }
//     }

//     write_line_tiff(output, mask_line, line);
//   }

//   for (int i = 0; i < 7; i++)
//   {
//     free(bufferTS[i]);
//     free(bufferNDVI[i]);
//     free(bufferAlb[i]);
//   }
//   free(bufferAlb);
//   free(bufferNDVI);
//   free(bufferTS);
// }

// /**
//  * @brief   Removes from the binary TIFF gave as input, groups of pixel with value 1, that have less pixels than a specified value.
//  * @param   input: A binary TIFF to be processed.
//  * @param   output: A binary TIFF.
//  * @param   groupSize: The inferior limit of the group of pixels size.
//  */
// void testMorphological(TIFF *input, TIFF *output, int groupSize)
// {

//   // Read the entire TIFF to the memory
//   // Create an same size matrix to serve as output

//   uint32 height_band, width_band;
//   TIFFGetField(input, TIFFTAG_IMAGELENGTH, &height_band);
//   TIFFGetField(input, TIFFTAG_IMAGEWIDTH, &width_band);

//   int **inputM = (int **)malloc(height_band * sizeof(int *));
//   int **outputM = (int **)malloc(height_band * sizeof(int *));

//   for (int i = 0; i < height_band; i++)
//   {

//     inputM[i] = (int *)malloc(width_band * sizeof(int));
//     outputM[i] = (int *)malloc(width_band * sizeof(int));
//   }

//   for (int i = 0; i < height_band; i++)
//   {

//     read_line_tiff(input, inputM[i], i);
//     read_line_tiff(input, outputM[i], i);
//   }

//   // Apply the routine

//   queue<pair<int, int>> fila;
//   set<pair<int, int>> cont;

//   for (int line = 0; line < height_band; line++)
//   {

//     for (int col = 0; col < width_band; col++)
//     {

//       if (inputM[line][col] == 1)
//       {

//         fila.push({line, col});
//         cont.insert({line, col});
//         inputM[line][col] = -1;

//         while (!fila.empty())
//         {

//           int i = fila.front().first;
//           int j = fila.front().second;
//           fila.pop();

//           if (j + 1 < width_band)
//           {

//             if (inputM[i][j + 1] == 1)
//             {
//               fila.push({i, j + 1});
//               cont.insert({i, j + 1});
//               inputM[i][j + 1] = -1;
//             }

//             if (i + 1 < height_band && inputM[i + 1][j + 1] == 1)
//             {
//               fila.push({i + 1, j + 1});
//               cont.insert({i + 1, j + 1});
//               inputM[i + 1][j + 1] = -1;
//             }

//             if (i > 0 && inputM[i - 1][j + 1] == 1)
//             {
//               fila.push({i - 1, j + 1});
//               cont.insert({i - 1, j + 1});
//               inputM[i - 1][j + 1] = -1;
//             }
//           }

//           if (j > 0)
//           {

//             if (inputM[i][j - 1] == 1)
//             {
//               fila.push({i, j - 1});
//               cont.insert({i, j - 1});
//               inputM[i][j - 1] = -1;
//             }

//             if (i + 1 < height_band && inputM[i + 1][j - 1] == 1)
//             {
//               fila.push({i + 1, j - 1});
//               cont.insert({i + 1, j - 1});
//               inputM[i + 1][j - 1] = -1;
//             }

//             if (i > 0 && inputM[i - 1][j - 1] == 1)
//             {
//               fila.push({i - 1, j - 1});
//               cont.insert({i - 1, j - 1});
//               inputM[i - 1][j - 1] = -1;
//             }
//           }

//           if (i + 1 < height_band && inputM[i + 1][j] == 1)
//           {
//             fila.push({i + 1, j});
//             cont.insert({i + 1, j});
//             inputM[i + 1][j] = -1;
//           }

//           if (i > 0 && inputM[i - 1][j] == 1)
//           {
//             fila.push({i - 1, j});
//             cont.insert({i - 1, j});
//             inputM[i - 1][j] = -1;
//           }
//         }

//         int group = cont.size();

//         for (auto elem : cont)
//         {

//           outputM[elem.first][elem.second] = (group >= groupSize);
//         }

//         cont.clear();
//       }
//       else if (inputM[line][col] == 0)
//       {

//         outputM[line][col] = 0;
//       }
//     }
//   }

//   // Write output TIFF

//   for (int i = 0; i < height_band; i++)
//   {

//     write_line_tiff(output, outputM[i], i);
//   }

//   for (int i = 0; i < height_band; i++)
//   {

//     free(inputM[i]);
//     free(outputM[i]);
//   }

//   free(inputM);
//   free(outputM);
// }

// /**
//  * @brief  Computes the HO.
//  * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
//  * @param  soil_heat_flux[]: Array containing the specified line from the G computation.
//  * @param  width_band: Band width.
//  * @param  ho_line[]: Auxiliary array for save the calculated value of HO for the line.
//  */
// void hoCalc(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[])
// {
//   for (int col = 0; col < width_band; col++)
//     ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];
// };