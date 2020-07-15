#include <complex.h>
#include <mex.h>

#include "../Setup/helper.h"

// Replace the calculation expV = exp((-1i*(dt/2)/hBar)*V);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long potential_pointer = mxGetScalar(prhs[0]);
	int nx = mxGetScalar(prhs[1]);
	int ny = mxGetScalar(prhs[2]);
	int nz = mxGetScalar(prhs[3]);
	double dt = mxGetScalar(prhs[4]);
	double h_bar = mxGetScalar(prhs[5]);

	// Pointer to a potential calculated in another function
	double *potential = int_to_pointer(potential_pointer);

	// Static pointer to the exponent
	static double *exp_v = NULL;
	double complex z1 = I * I;
	if (exp_v == NULL) {
		cudaMallocManaged(&exp_v, nx * ny * nz * sizeof(bla));
	}

	// Create output variables to store the pointers to the potentials and avoid memory leaks
	plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	long long *ptr_exp_v = (long long *) mxGetPr(plhs[0]);

	// Set the union to the value of the pointer and the return array to the correesponding
	// integer value.
	ptr_exp_v[0] = pointer_to_int(exp_v);
}
