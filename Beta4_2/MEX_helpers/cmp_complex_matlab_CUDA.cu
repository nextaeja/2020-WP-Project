#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <stdlib.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if (!mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("Wavepacket:ArrayCmp", "MATLAB array must be complex\n");
	}
	double *matlab_array_real = mxGetPr(prhs[0]);
	double *matlab_array_imag = mxGetPi(prhs[0]);
	long long cuda_array_ptr = mxGetScalar(prhs[1]);
	double epsilon = mxGetScalar(prhs[2]);
	int nx = mxGetScalar(prhs[3]);
	int ny = mxGetScalar(prhs[4]);
	int nz = mxGetScalar(prhs[5]);
	int debug = 0;
	if (nrhs > 6) {
		debug = mxGetScalar(prhs[6]);
	}
	size_t size = nx*ny*nz;

	cuComplex *dev_cuda_array = reinterpret_cast<cuComplex *>(cuda_array_ptr);
	cuComplex *cuda_array = reinterpret_cast<cuComplex *>(malloc(size * sizeof(cuComplex)));
	cudaMemcpy(cuda_array, dev_cuda_array, size * sizeof(cuComplex), cudaMemcpyDeviceToHost);

	for (int k=0; k<nz; k++) {
		for (int j=0; j<ny; j++) {
			for (int i=0; i<nx; i++) {
				int idx = i + nx*j + nx*ny*k;
				double real_diff = abs((matlab_array_real[idx] - cuda_array[idx].x) / matlab_array_real[idx]);
				double imag_diff = abs((matlab_array_imag[idx] - cuda_array[idx].y) / matlab_array_imag[idx]);

				if (debug) {
					mexPrintf("(%d %d %d %d) (%e %e) (%e %e) (%e %e)\n", i, j, k, idx, matlab_array_real[idx], matlab_array_imag[idx], cuda_array[idx].x, cuda_array[idx].y, real_diff, imag_diff);
				}

				if (real_diff > epsilon || imag_diff > epsilon) {
					mexPrintf("(%d %d %d %d) (%e %e) (%e %e) (%e %e)\n", i, j, k, idx, matlab_array_real[idx], matlab_array_imag[idx], cuda_array[idx].x, cuda_array[idx].y, real_diff, imag_diff);
					mexErrMsgIdAndTxt("Wavepacket:ArrayCmp", "Arrays are not equal\n");
				}
			}
		}
	}

	free(cuda_array);
}
