/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <assert.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

#define N_POINTERS 8

__global__ void initialize_array(double *array, size_t size) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < size) {
        array[tid] = 0.0;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void initialize_array(myComplex *array, size_t size) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < size) {
        myComplex w;
        w.x = 0.0;
        w.y = 0.0;
        array[tid] = w;

        tid += blockDim.x * gridDim.x;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	size_t nx = mxGetScalar(prhs[0]);
	size_t ny = mxGetScalar(prhs[1]);
	size_t nz = mxGetScalar(prhs[2]);
	int num_adsorbates = mxGetScalar(prhs[3]);
    int num_iterations = mxGetScalar(prhs[4]);
	num_adsorbates = (num_adsorbates < 1) ? 1 : num_adsorbates;

	double *z_offset;
	double *dev_x0;
	double *dev_y0;
	double *kSquared;
    double *gaussian_positions;
	myComplex *exp_v;
	myComplex *exp_k;
    myComplex *psi;

	// This data is real
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&z_offset), nx * ny * sizeof(double)));
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&dev_x0), num_adsorbates * sizeof(double)));
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&dev_y0), num_adsorbates * sizeof(double)));
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&kSquared), nx * ny * nz * sizeof(double)));
    CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&gaussian_positions), num_adsorbates * 2 * num_iterations * sizeof(double)));

	// This data is complex
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&exp_v), nx * ny * nz * sizeof(myComplex)));
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&exp_k), nx * ny * nz * sizeof(myComplex)));
    CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&psi), nx * ny * nz * sizeof(myComplex)));

	// Initialize all arrays
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(z_offset, nx * ny);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(dev_x0, num_adsorbates);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(dev_y0, num_adsorbates);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(kSquared, nx * ny * nz);
    initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(gaussian_positions, num_adsorbates * 2 * num_iterations);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(exp_v, nx * ny * nz);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(exp_k, nx * ny * nz);
    initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(psi, nx * ny * nz);

	plhs[0] = mxCreateNumericMatrix(1, N_POINTERS, mxINT64_CLASS, mxREAL);
	long long *ptr_potential = (long long *) mxGetPr(plhs[0]);

	ptr_potential[0] = reinterpret_cast<long long>(z_offset);
	ptr_potential[1] = reinterpret_cast<long long>(dev_x0);
	ptr_potential[2] = reinterpret_cast<long long>(dev_y0);
	ptr_potential[3] = reinterpret_cast<long long>(kSquared);
	ptr_potential[4] = reinterpret_cast<long long>(exp_v);
	ptr_potential[5] = reinterpret_cast<long long>(exp_k);
    ptr_potential[6] = reinterpret_cast<long long>(psi);
    ptr_potential[7] = reinterpret_cast<long long>(gaussian_positions);
}
