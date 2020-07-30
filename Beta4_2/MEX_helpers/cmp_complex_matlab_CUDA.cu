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
	int size = mxGetScalar(prhs[3]);
	int debug = 0;
	if (nrhs > 4) {
		debug = mxGetScalar(prhs[4]);
	}

	cuComplex *dev_cuda_array = reinterpret_cast<cuComplex *>(cuda_array_ptr);
	cuComplex *cuda_array = reinterpret_cast<cuComplex *>(malloc(size * sizeof(cuComplex)));
	cudaMemcpy(cuda_array, dev_cuda_array, size * sizeof(cuComplex), cudaMemcpyDeviceToHost);

	for (int i=0; i<size; i++) {
		double real_diff = abs(matlab_array_real[i] - cuda_array[i].x);
		double imag_diff = abs(matlab_array_imag[i] - cuda_array[i].y);

		if (debug) {
			mexPrintf("(%e %e) (%e %e) (%e %e)\n", matlab_array_real[i], matlab_array_imag[i], cuda_array[i].x, cuda_array[i].y, real_diff, imag_diff);
		}

		if (real_diff > epsilon || imag_diff > epsilon) {
			mexErrMsgIdAndTxt("Wavepacket:ArrayCmp", "Arrays are not equal\n");
		}
	}

	free(cuda_array);
}
