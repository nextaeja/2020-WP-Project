#include <mex.h>
#include <matrix.h>
#include <math.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"

__global__ void copy_complex_array(complex *dest, double *real, double *imag, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		dest[tid] = complex(real[tid], imag[tid]);

		tid += blockDim.x * gridDim.x;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	long long dest_ptr = mxGetScalar(prhs[0]);
	double *source_real = mxGetPr(prhs[1]);
	double *source_imag = mxGetPi(prhs[1]);
	size_t size = mxGetScalar(prhs[2]);

	complex *dest = reinterpret_cast<complex *>(dest_ptr);

	// Allocate the space on the GPU
	double *dev_source_real, *dev_source_imag;
	cudaMallocManaged(reinterpret_cast<void **>(&dev_source_real), size * sizeof(double));
	cudaMallocManaged(reinterpret_cast<void **>(&dev_source_imag), size * sizeof(double));

	// Copy input data to GPU
	cudaMemcpy(dev_source_real, source_real, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_source_imag, source_imag, size * sizeof(double), cudaMemcpyHostToDevice);

	copy_complex_array<<<NUM_BLOCKS, NUM_THREADS>>>(dest, dev_source_real, dev_source_imag, size);

	cudaFree(dev_source_real);
	cudaFree(dev_source_imag);
}
