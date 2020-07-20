#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"

const int NUM_BLOCKS = 65536;
const int NUM_THREADS = 256;

// Compute (-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass))
__global__ void shift_ksquared(complex *expk, double *ksquared, complex prefactor, double h_bar, double mass, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i=tid; i<size; i += num_threads) {
		expk[i] = prefactor * (h_bar*h_bar * ksquared[i] / (2*mass));
	}
}

__global__ void bangalla()

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

	// Copy k_squared over to GPU memory
	cudaMemcpy(reinterpret_cast<void **>(&dev_k_squared),
			   reinterpret_cast<void **>(&k_squared),
		   	   size * sizeof(double),
		       cudaMemcpyHostToDevice);

	// Get scaling constant (-1i*dt/hBar)
	complex scale = -complex(0,1)*(dt / h_bar);

	// Compute the argument of the exponential
	//shift_ksquared<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expk, dev_k_squared, scale, h_bar, mass, size);

	// Compute the exponential
	//cexp<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expk, size);
}
