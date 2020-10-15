/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

__global__ void copy_complex_array(myComplex *source, double *real, double *imag, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		real[tid] = source[tid].x;
		imag[tid] = source[tid].y;

		tid += blockDim.x * gridDim.x;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long source_ptr = mxGetScalar(prhs[0]);
	size_t nx = mxGetScalar(prhs[1]);
	size_t ny = mxGetScalar(prhs[2]);
	size_t nz = mxGetScalar(prhs[3]);
	size_t size = nx * ny * nz;

	myComplex *dev_source = reinterpret_cast<myComplex *>(source_ptr);

	// Setup auxiliary variables
	double *dev_dest_real, *dev_dest_imag;
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&dev_dest_real), size * sizeof(double)));
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&dev_dest_imag), size * sizeof(double)));

	// Copy from myComplex to arrays of double
	copy_complex_array<<<NUM_BLOCKS, NUM_THREADS>>>(dev_source, dev_dest_real, dev_dest_imag, size);

	// Set the output matrix
	const mwSize dims[] = {nx, ny, nz};
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxCOMPLEX);
	double *dest_real = mxGetPr(plhs[0]);
	double *dest_imag = mxGetPi(plhs[0]);

	// Copy from device to host
	CUDA_HANDLE(cudaMemcpy(dest_real, dev_dest_real, size * sizeof(double), cudaMemcpyDeviceToHost));
	CUDA_HANDLE(cudaMemcpy(dest_imag, dev_dest_imag, size * sizeof(double), cudaMemcpyDeviceToHost));

	// Free auxiliary
	CUDA_HANDLE(cudaFree(dev_dest_real));
	CUDA_HANDLE(cudaFree(dev_dest_imag));
}
