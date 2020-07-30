#ifndef __COMPLEX__
#define __COMPLEX__

#include <math.h>

typedef float2 cuComplex;

// Overloaded addition
__device__ cuComplex element_complex_add(cuComplex z1, cuComplex z2);
__device__ cuComplex element_complex_add(float x, cuComplex z);
__device__ cuComplex element_complex_add(cuComplex z, float x);

// Overloaded multiplication
__device__ cuComplex element_complex_mul(cuComplex z1, cuComplex z2);
__device__ cuComplex element_complex_mul(float x, cuComplex z);
__device__ cuComplex element_complex_mul(cuComplex z, float x);

__device__ cuComplex element_complex_exp(cuComplex z);

__global__ void complex_scale_shift(cuComplex *zs, float scale, float shift, bool opposite, size_t size);
__global__ void complex_exp(cuComplex *zs, size_t size);
__global__ void complex_mul(cuComplex *zs, cuComplex *ws, size_t size);

#endif
