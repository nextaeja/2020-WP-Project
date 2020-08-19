#include <stdio.h>

#include "complex.h"

/// =============== SCALAR FUNCTIONS ===============
__device__ myComplex element_complex_add(myComplex z1, myComplex z2) {
	myComplex w;
	w.x = z1.x + z2.x;
	w.y = z1.y + z2.y;

	return w;
}

__device__ myComplex element_complex_add(double x, myComplex z) {
	myComplex w;
	w.x = x + z.x;
	w.y = z.y;

	return w;
}

__device__ myComplex element_complex_add(myComplex z, double x) {
	return element_complex_add(x, z);
}

__device__ myComplex element_complex_mul(myComplex z1, myComplex z2) {
	myComplex w;
	w.x = z1.x*z2.x - z1.y*z2.y;
	w.y = z1.x*z2.y + z1.y*z2.x;

	return w;
}

__device__ myComplex element_complex_mul(double x, myComplex z) {
	myComplex w;
	w.x = x * z.x;
	w.y = x * z.y;

	return w;
}

__device__ myComplex element_complex_mul(myComplex z, double x) {
	return element_complex_mul(x, z);
}

__device__ myComplex element_complex_exp(myComplex z) {
	myComplex w;
	w.x = cos(z.y);
	w.x = sin(z.y);

	return element_complex_mul(exp(z.x), w);
}
/// ================================================

/// =============== VECTOR FUNCTIONS ===============
// Scale and shift an array by real constants
// i.e.: z --> scale*z - shift	if opposite is FALSE
//	     z --> shift - scale*z  if opposite is TRUE
__global__ void complex_scale_shift(myComplex *zs, double scale, double shift, bool opposite, size_t size) {
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

// Scale an array by real constants
// i.e.: z --> scale*z - shift	if opposite is FALSE
__global__ void complex_scale(myComplex *zs, double scale, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		zs[tid] = element_complex_mul(scale, zs[tid]);

		tid += blockDim.x * gridDim.x;
	}
}

// Compute the exponential of an array
// 		z --> exp(z)
__global__ void complex_exp(myComplex *zs, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		zs[tid] = element_complex_exp(zs[tid]);

		tid += blockDim.x * gridDim.x;
	}
}

// Multiply two arrays elementwise
//		(x, y) --> x*y
__global__ void complex_mul(myComplex *zs, myComplex *ws, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		zs[tid] = element_complex_mul(zs[tid], ws[tid]);

		tid += blockDim.x * gridDim.x;
	}
}
/// ================================================
