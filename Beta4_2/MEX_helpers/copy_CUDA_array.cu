/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/cuda_helper.h"

__global__ void copy_array(double *dest, double *source, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		dest[tid] = source[tid];

		tid += blockDim.x * gridDim.x;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long dest_ptr = mxGetScalar(prhs[0]);
	double *source = mxGetPr(prhs[1]);
	size_t size = mxGetScalar(prhs[2]);

	double *dest = reinterpret_cast<double *>(dest_ptr);

	// Allocate the space on the GPU
	double *dev_source;
	cudaMallocManaged(reinterpret_cast<void **>(&dev_source), size * sizeof(double));

	// Copy input data to GPU
	cudaMemcpy(dev_source, source, size * sizeof(double), cudaMemcpyHostToDevice);

	copy_array<<<NUM_BLOCKS, NUM_THREADS>>>(dest, dev_source, size);

	cudaFree(dev_source);
}
