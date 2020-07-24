#ifndef __COMPLEX__
#define __COMPLEX__

#include <math.h>

typedef float2 Complex;

// Overloaded addition
static __device__ __host__ inline Complex _complex_add(Complex z1, Complex z2);
static __device__ __host__ inline Complex _complex_add(float x, Complex z);
static __device__ __host__ inline Complex _complex_add(Complex z, float x);

// Overloaded multiplication
static __device__ __host__ inline Complex _complex_mul(Complex z1, Complex z2);
static __device__ __host__ inline Complex _complex_mul(float x, Complex z);
static __device__ __host__ inline Complex _complex_mul(Complex z, float x);

static __device__ __host__ inline Complex _complex_exp(Complex z);

__global__ void complex_scale_shift(Complex *zs, float scale, float shift, bool opposite, size_t size);
__global__ void complex_exp(Complex *zs, size_t size);

#endif
