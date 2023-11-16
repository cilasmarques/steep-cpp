#include "utils.h"

__global__ void correctionCycle(double *surfaceTemperatureLine, double *zomLine, double *ustarRLine, double *ustarWLine, double *rahRLine,
                                double *rahWLine, double *a, double *b, double *u200, int *size);

__global__ void correctionCycleSTEEP(double *surfaceTemperatureLine, double *d0Line, double *kb1_vector, double *zomLine, double *ustarRLine,
                                     double *ustarWLine, double *rahRLine, double *rahWLine, double *a, double *b, double *u200, int *size);