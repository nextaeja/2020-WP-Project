#include <stdio.h>

#include "complex.h"

/// =============== SCALAR FUNCTIONS ===============
__device__ cuComplex element_complex_add(cuComplex z1, cuComplex z2) {
	cuComplex w;
	w.x = z1.x + z2.x;
	w.y = z1.y + z2.y;

	return w;
}

__device__ cuComplex element_complex_add(float x, cuComplex z) {
	cuComplex w;
	w.x = x + z.x;
	w.y = z.y;

	return w;
}

__device__ cuComplex element_complex_add(cuComplex z, float x) {
	return element_complex_add(x, z);
}

__device__ cuComplex element_complex_mul(cuComplex z1, cuComplex z2) {
	cuComplex w;
	w.x = z1.x*z2.x - z1.y*z2.y;
	w.y = z1.x*z2.y + z1.y*z2.x;

	return w;
}

__device__ cuComplex element_complex_mul(float x, cuComplex z) {
	cuComplex w;
	w.x = x * z.x;
	w.y = x * z.y;

	return w;
}

__device__ cuComplex element_complex_mul(cuComplex z, float x) {
	return element_complex_mul(x, z);
}

__device__ cuComplex element_complex_exp(cuComplex z) {
	cuComplex w;
	w.x = cos(z.y);
	w.x = sin(z.y);

	return element_complex_mul(exp(z.x), w);
}
/// ================================================

/// =============== VECTOR FUNCTIONS ===============
// Scale and shift an array by real constants
// i.e.: z --> scale*z - shift	if opposite is FALSE
//	     z --> shift - scale*z  if opposite is TRUE
__global__ void complex_scale_shift(cuComplex *zs, float scale, float shift, bool opposite, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		if (opposite) {
			zs[tid] = element_complex_add(element_complex_mul(scale, zs[tid]), -shift);
		} else {
			zs[tid] = element_complex_add(shift, element_complex_mul(-scale, zs[tid]));
		}

		tid += blockDim.x * gridDim.x;
	}
}

// Compute the exponential of an array
// 		z --> exp(z)
__global__ void complex_exp(cuComplex *zs, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		zs[tid] = element_complex_exp(zs[tid]);

		tid += blockDim.x * gridDim.x;
	}
}

// Multiply two arrays elementwise
//		(x, y) --> x*y
__global__ void complex_mul(cuComplex *zs, cuComplex *ws, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		zs[tid] = element_complex_mul(zs[tid], ws[tid]);

		tid += blockDim.x * gridDim.x;
	}
}
/// ================================================
