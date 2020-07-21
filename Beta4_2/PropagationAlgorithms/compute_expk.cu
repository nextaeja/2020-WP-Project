#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

// Compute (-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass))
__global__ void shift_ksquared(complex *expk, double *ksquared, complex prefactor, double h_bar, double mass, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		expk[tid] = prefactor * (h_bar*h_bar * ksquared[tid] / (2*mass));

		tid += blockDim.x * gridDim.x;
	}
}

// Compute exp((-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass)))
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Exctract the parameters
	long long expk_ptr = mxGetScalar(prhs[0]);
	long long k_squared_ptr = mxGetScalar(prhs[1]);
	double h_bar = mxGetScalar(prhs[2]);
	double dt = mxGetScalar(prhs[3]);
	double mass = mxGetScalar(prhs[4]);
	double *k_squared = mxGetPr(prhs[5]);
	size_t size = mxGetScalar(prhs[6]);

	// Parse the pointer to allocated space for expK and k_squared
	complex *dev_expk = reinterpret_cast<complex *>(expk_ptr);
	double *dev_k_squared = reinterpret_cast<double *>(k_squared_ptr);

	// k_squared is computed in MATLAB and copied over the GPU
	cudaMemcpy(dev_k_squared, k_squared, size * sizeof(double), cudaMemcpyHostToDevice);

	// Get scaling constant (-1i*dt/hBar)
	complex scale = -complex(0,1)*(dt / h_bar);

	// Compute the argument of the exponential
	shift_ksquared<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expk, dev_k_squared, scale, h_bar, mass, size);

	// Compute the exponential
	cexp<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expk, size);
}
