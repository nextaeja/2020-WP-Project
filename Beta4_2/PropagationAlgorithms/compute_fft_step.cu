#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

__global__ void psi_v_half_step(Complex *psi, Complex *expv, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		psi[tid] = _complex_mul(psi[tid], expv[tid]);

		tid += blockDim.x * gridDim.x;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long expv_ptr = mxGetScalar(prhs[0]);
	long long expk_ptr = mxGetScalar(prhs[1]);
	long long psi_ptr = mxGetScalar(prhs[2]);
	int nx = mxGetScalar(prhs[3]);
	int ny = mxGetScalar(prhs[4]);
	int nz = mxGetScalar(prhs[5]);

	// Parse the pointers
	Complex *dev_expv = reinterpret_cast<Complex *>(expv_ptr);
	Complex *dev_expk = reinterpret_cast<Complex *>(expk_ptr);
	Complex *dev_psi = reinterpret_cast<Complex *>(psi_ptr);
}
