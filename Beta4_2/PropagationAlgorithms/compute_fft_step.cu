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
	myComplex *dev_expv = reinterpret_cast<myComplex *>(expv_ptr);
	myComplex *dev_expk = reinterpret_cast<myComplex *>(expk_ptr);
	myComplex *dev_psi = reinterpret_cast<myComplex *>(psi_ptr);

	// Plan the FFT
	cufftHandle forward_plan, inverse_plan;
	int n[3] = {nz, ny, nx};
	int idist = size;
	int odist = size;
	int istride = 1;
	int ostride = 1;
	int inembed[3] = {nz, ny, nx}; // MATLAB inverts rows and columns
	int onembed[3] = {nz, ny, nx};
	cufftPlanMany(&forward_plan, 3, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, 1);
	cufftPlanMany(&inverse_plan, 3, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2Z, 1);

	// psiVStepHalf = expV.*psi;
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expv, size);
	cudaDeviceSynchronize();

	// Compute the forward FFT
	cufftExecZ2Z(forward_plan, dev_psi, dev_psi, CUFFT_FORWARD);

	// psiKStepFT = expK.*psiVStepHalfFT;
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expk, size);

	// Invert FFT
	cufftExecZ2Z(inverse_plan, dev_psi, dev_psi, CUFFT_INVERSE);
	complex_scale<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, 1/(double) size, size);

	cufftDestroy(forward_plan);
	cufftDestroy(inverse_plan);
}
