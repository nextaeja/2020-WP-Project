#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

#define NDIMS 3

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
	int n[3] = {nx, ny, nz};
	cufftPlan3d(&plan, nx, ny, nz, CUFFT_C2C);
	/*
	cufftPlanMany(&plan, 3, n,
					NULL, 1, size,
					NULL, 1, size,
					CUFFT_C2C, 1);
	*/
	cufftExecC2C(plan, dev_psi, dev_psi, CUFFT_INVERSE);

	/*
	// psiKStepFT = expK.*psiVStepHalfFT;
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expk, size);

	// Invert FFT
	cufftExecC2C(plan, dev_psi, dev_psi, CUFFT_INVERSE);
	*/

	cufftDestroy(plan);
}
