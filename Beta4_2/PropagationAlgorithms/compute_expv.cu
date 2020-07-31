#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

// Compute exp((-1i*(dt/2)/hBar)*V)
__global__ void compute_expv(myComplex *potential, double scale, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		myComplex expv;
		expv.x = cos(scale * potential[tid].x);
		expv.y = sin(scale * potential[tid].x);

		potential[tid] = expv;

		tid += blockDim.x * gridDim.x;
	}
}

// expV = exp((-1i*(dt/2)/hBar)*V);
// Compute the exponential of the potential in place (no copy/allocation is needed)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Exctract the parameters
	long long potential_ptr = mxGetScalar(prhs[0]);
	double h_bar = mxGetScalar(prhs[1]);
	double dt = mxGetScalar(prhs[2]);
	size_t size = mxGetScalar(prhs[3]);

	// Parse the pointer to allocated space for expK and k_squared
	myComplex *potential = reinterpret_cast<myComplex *>(potential_ptr);

	// Get scaling constant (-1i*(dt/2)/hBar)
	double scale = -(dt/2) / h_bar;

	// Calculate the exponential
	compute_expv<<<NUM_BLOCKS, NUM_THREADS>>>(potential, scale, size);
}
