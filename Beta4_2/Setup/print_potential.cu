#include <mex.h>

#include "helper.h"

// Print a 3d potential
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long potential_pointer = mxGetScalar(prhs[0]);
	int nx = mxGetScalar(prhs[1]);
	int ny = mxGetScalar(prhs[2]);
	int nz = mxGetScalar(prhs[3]);

	double *potential = int_to_pointer(potential_pointer);

	for (int k=0; k<nz; k++) {
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				int idx = k*nx*ny+j*nx+i;
				mexPrintf("%e ", potential[idx]);
			}
			mexPrintf("\n");
		}
		mexPrintf("\n");
	}
}
