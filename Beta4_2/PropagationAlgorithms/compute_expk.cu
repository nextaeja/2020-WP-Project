/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

// Compute (-1i*dt/hBar)*(hBar^2*kSquared/(2*mass))
__global__ void compute_expk(myComplex *expk, double *ksquared, double dt, double h_bar, double mass, size_t size) {
	myComplex element;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		double argument = h_bar*h_bar * ksquared[tid] / (2*mass);
		argument *= -dt / h_bar;

		element.x = cos(argument);
		element.y = sin(argument);

		expk[tid] = element;

		tid += blockDim.x * gridDim.x;
	}
}

// Compute exp((-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass)))
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Exctract the parameters
	long long expk_ptr = mxGetScalar(prhs[0]);
	long long k_squared_ptr = mxGetScalar(prhs[1]);
	double h_bar = mxGetScalar(prhs[2]);
	double dt = mxGetScalar(prhs[3]);
	double mass = mxGetScalar(prhs[4]);
	double *k_squared = mxGetPr(prhs[5]);
	size_t size = mxGetScalar(prhs[6]);

	// Parse the pointer to allocated space for expK and k_squared
	myComplex *dev_expk = reinterpret_cast<myComplex *>(expk_ptr);
	double *dev_k_squared = reinterpret_cast<double *>(k_squared_ptr);

	// k_squared is computed in MATLAB and copied over the GPU
	CUDA_HANDLE(cudaMemcpy(dev_k_squared, k_squared, size * sizeof(double), cudaMemcpyHostToDevice));

	// Compute the argument of the exponential
	compute_expk<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expk, dev_k_squared, dt, h_bar, mass, size);
}
