#include <stdio.h>

#include "complex.h"

// Exponential of an array of complex numbers
__global__ void cexp(complex *zs, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i=tid; i<size; i += num_threads) {
		zs[i] = zs[i]._cexp();
	}
}
