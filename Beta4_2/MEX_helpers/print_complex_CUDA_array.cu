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

	Complex *dev_potential = reinterpret_cast<Complex *>(potential_ptr);
	Complex *potential = reinterpret_cast<Complex *>(malloc(nx * ny * nz * sizeof(complex)));
	cudaMemcpy(potential, dev_potential, nx * ny * nz * sizeof(complex), cudaMemcpyDeviceToHost);

	for (int k=0; k<nz; k++) {
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				int idx = k*nx*ny+j*nx+i;
				mexPrintf("(%e + i*%e) ", potential[idx].x, potential[idx].y);
			}
			mexPrintf("\n");
		}
		mexPrintf("\n");
	}

	free(potential);
}
