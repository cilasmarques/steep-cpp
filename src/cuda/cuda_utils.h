#pragma once

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "cuda_profiler_api.h"
#include <cuda_runtime_api.h>

/**
 * This function checks the return value of the CUDA runtime call and exits
 * the application if the call failed.
 *
 * @param err: CUDA error code.
 * @param file: File name where the error occurred.
 * @param line: Line number where the error occurred.
 */
static void HandleError(cudaError_t err, const char *file, int line)
{
  if (err != cudaSuccess)
  {
    printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
}

/**
 * This function checks the return value of the CUDA runtime call and exits
 * the application if the call failed.
 *
 * @param err: CUDA error code.
 * @param file: File name where the error occurred.
 * @param line: Line number where the error occurred.
 */
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

/**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 *
 * See cuda.h for error code descriptions.
 */
#define CHECK_CUDA_RESULT(N)                                       \
  {                                                                \
    CUresult result = N;                                           \
    if (result != 0)                                               \
    {                                                              \
      printf("CUDA call on line %d returned error %d\n", __LINE__, \
             result);                                              \
      exit(1);                                                     \
    }                                                              \
  }
