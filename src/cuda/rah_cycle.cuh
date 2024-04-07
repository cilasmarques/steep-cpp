#include "utils.h"

__global__ void correctionCycle(double *surfaceTemperatureLine, double *zomLine, double *ustarRLine, double *ustarWLine, double *rahRLine,
                                double *rahWLine, double *a, double *b, double *u200, int *size);

__global__ void correctionCycleSTEEP(double *surface_temperature_pointer, double *d0_pointer, double *kb1_pointer, double *zom_pointer, double *ustarR_pointer,
                                     double *ustarW_pointer, double *rahR_pointer, double *rahWL_pointer, double *H_pointer, double a, double b, int height,
                                     int width);
