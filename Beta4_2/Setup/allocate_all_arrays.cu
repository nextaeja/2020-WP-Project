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

__global__ void initialize_array(complex *array, size_t size) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < size) {
        array[tid] = complex(0.0, 0.0);

        tid += blockDim.x * gridDim.x;
    }
}

/*
 * Allocate memory for:
 *	- potential (cuda_setup_dynamic_potential_2.cu)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	size_t nx = mxGetScalar(prhs[0]);
	size_t ny = mxGetScalar(prhs[1]);
	size_t nz = mxGetScalar(prhs[2]);
	int num_adsorbates = mxGetScalar(prhs[3]);
	num_adsorbates = (num_adsorbates < 1) ? 1 : num_adsorbates;

	double *z_offset;
	double *dev_x0;
	double *dev_y0;
	double *kSquared;
    complex *potential;
	complex *exp_v;
	complex *exp_k;
    complex *psi_ptr;

	// This data is real
	cudaMallocManaged(reinterpret_cast<void **>(&z_offset), nx * ny * sizeof(double));
	cudaMallocManaged(reinterpret_cast<void **>(&dev_x0), num_adsorbates * sizeof(double));
	cudaMallocManaged(reinterpret_cast<void **>(&dev_y0), num_adsorbates * sizeof(double));
	cudaMallocManaged(reinterpret_cast<void **>(&kSquared), nx * ny * nz * sizeof(double));

	// This data is complex
    cudaMallocManaged(reinterpret_cast<void **>(&potential), nx * ny * nz * sizeof(complex));
	cudaMallocManaged(reinterpret_cast<void **>(&exp_v), nx * ny * nz * sizeof(complex));
	cudaMallocManaged(reinterpret_cast<void **>(&exp_k), nx * ny * nz * sizeof(complex));
    cudaMallocManaged(reinterpret_cast<void **>(&psi_ptr), nx * ny * nz * sizeof(complex));

	// Initialize all arrays
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(z_offset, nx * ny);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(dev_x0, num_adsorbates);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(dev_y0, num_adsorbates);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(kSquared, nx * ny * nz);
    initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(potential, nx * ny * nz);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(exp_v, nx * ny * nz);
	initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(exp_k, nx * ny * nz);
    initialize_array<<<NUM_BLOCKS, NUM_THREADS>>>(psi_ptr, nx * ny * nz);

	plhs[0] = mxCreateNumericMatrix(1, N_POINTERS, mxINT64_CLASS, mxREAL);
	long long *ptr_potential = (long long *) mxGetPr(plhs[0]);

	ptr_potential[0] = reinterpret_cast<long long>(potential);
	ptr_potential[1] = reinterpret_cast<long long>(z_offset);
	ptr_potential[2] = reinterpret_cast<long long>(dev_x0);
	ptr_potential[3] = reinterpret_cast<long long>(dev_y0);
	ptr_potential[4] = reinterpret_cast<long long>(kSquared);
	ptr_potential[5] = reinterpret_cast<long long>(exp_v);
	ptr_potential[6] = reinterpret_cast<long long>(exp_k);
    ptr_potential[7] = reinterpret_cast<long long>(psi_ptr);
}
