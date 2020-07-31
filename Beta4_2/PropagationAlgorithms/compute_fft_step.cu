#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

#define NDIMS 3

__global__ void swap_axes(cuComplex *dev_psi, int nx, int ny, int nz) {
	for (int k=0; k<nz; k++) {
		for (int j=0; j<ny; j++) {
			for (int i=0; i<nx; i++) {
				int source_idx = i + j*nx + k*nx*ny;
				int dest_idx = j + i*ny + k*nx*ny;

				cuComplex tmp = dev_psi[dest_idx];
				dev_psi[dest_idx] = dev_psi[source_idx];
				dev_psi[source_idx] = tmp;
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long expv_ptr = mxGetScalar(prhs[0]);
	long long expk_ptr = mxGetScalar(prhs[1]);
	long long psi_ptr = mxGetScalar(prhs[2]);
	int nx = mxGetScalar(prhs[3]);
	int ny = mxGetScalar(prhs[4]);
	int nz = mxGetScalar(prhs[5]);
	size_t size = nx * ny * nz;

	// Parse the pointers
	cuComplex *dev_expv = reinterpret_cast<cuComplex *>(expv_ptr);
	cuComplex *dev_expk = reinterpret_cast<cuComplex *>(expk_ptr);
	cuComplex *dev_psi = reinterpret_cast<cuComplex *>(psi_ptr);

	// psiVStepHalf = expV.*psi;
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expv, size);
	cudaDeviceSynchronize();

	// Compute the FFT
	cufftHandle plan;
	int n[3] = {nz, nx, ny};
	int idist = size;
	int odist = size;
	int istride = 1;
	int ostride = 1;
	int inembed[3] = {nz, nx, ny};
	int onembed[3] = {nz, nx, ny};
	cufftPlanMany(
		&plan, 3, n,
		inembed, istride, idist,
		onembed, ostride, odist,
		CUFFT_C2C, 1
	);
	cufftExecC2C(plan, dev_psi, dev_psi, CUFFT_FORWARD);

	/*
	// psiKStepFT = expK.*psiVStepHalfFT;
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expk, size);

	// Invert FFT
	cufftExecC2C(plan, dev_psi, dev_psi, CUFFT_INVERSE);
	*/

	cufftDestroy(plan);
}
