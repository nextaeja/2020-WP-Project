#include <stdio.h>

#include "complex.h"

/// =============== SCALAR FUNCTIONS ===============
__device__ __host__ inline Complex _complex_add(Complex z1, Complex z2) {
	Complex w;
	w.x = z1.x + z2.x;
	w.y = z1.y + z2.y;

	return w;
}

__device__ __host__ inline Complex _complex_add(float x, Complex z) {
	Complex w;
	w.x = x + z.x;
	w.y = z.y;

	return w;
}

__device__ __host__ inline Complex _complex_add(Complex z, float x) {
	return _complex_add(x, z);
}

__device__ __host__ inline Complex _complex_mul(Complex z1, Complex z2) {
	Complex w;
	w.x = z1.x*z2.x - z1.y*z2.y;
	w.y = z1.x*z2.y + z1.y*z2.x;

	return w;
}

__device__ __host__ inline Complex _complex_mul(float x, Complex z) {
	Complex w;
	w.x = x * z.x;
	w.y = x * z.y;

	return w;
}

__device__ __host__ inline Complex _complex_mul(Complex z, float x) {
	return _complex_mul(x, z);
}

__device__ __host__ inline Complex _complex_exp(Complex z) {
	Complex w;
	w.x = cos(z.y);
	w.x = sin(z.y);

	return _complex_mul(exp(z.x), w);
}
/// ================================================

/// =============== VECTOR FUNCTIONS ===============
// Scale and shift an array by real constants
// i.e.: z --> scale*z - shift	if opposite is FALSE
//	     z --> shift - scale*z  if opposite is TRUE
__global__ void complex_scale_shift(Complex *zs, float scale, float shift, bool opposite, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		if (opposite) {
			zs[tid] = _complex_add(_complex_mul(scale, zs[tid]), -shift);
		} else {
			zs[tid] = _complex_add(shift, _complex_mul(-scale, zs[tid]));
		}

		tid += blockDim.x * gridDim.x;
	}
}

// Compute the exponential of an array
// 		z --> exp(z)
__global__ void complex_exp(Complex *zs, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		zs[tid] = _complex_exp(zs[tid]);

		tid += blockDim.x * gridDim.x;
	}
}
/// ================================================
