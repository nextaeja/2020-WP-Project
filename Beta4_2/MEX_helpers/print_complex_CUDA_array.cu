#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long potential_ptr = mxGetScalar(prhs[0]);
	int nx = mxGetScalar(prhs[1]);
	int ny = mxGetScalar(prhs[2]);
	int nz = mxGetScalar(prhs[3]);

	complex *dev_potential = reinterpret_cast<complex *>(potential_ptr);
	complex *potential = reinterpret_cast<complex *>(malloc(nx * ny * nz * sizeof(complex)));
	cudaMemcpy(potential, dev_potential, nx * ny * nz * sizeof(complex), cudaMemcpyDeviceToHost);

	for (int k=0; k<nz; k++) {
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				int idx = k*nx*ny+j*nx+i;
				mexPrintf("(%e + i*%e) ", potential[idx].real, potential[idx].imag);
			}
			mexPrintf("\n");
		}
		mexPrintf("\n");
	}

	free(potential);
}
