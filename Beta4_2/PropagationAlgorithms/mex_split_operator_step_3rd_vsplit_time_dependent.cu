#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"
#include "../MEX_helpers/interpolation1d.h"

#define NDIMS 3


__global__ void compute_expv(myComplex *dev_expv, double scale, size_t size);

void split_operator_3rd_vsplit_time(myComplex *dev_psi, myComplex *dev_expv, myComplex *dev_expk, double *dev_gauss_time,
		double *dev_gauss_pos, double *dev_x0, double *dev_y0, double t_query, double A, double eV, double expv_scale, size_t size,
		cufftHandle forward_plan, cufftHandle inverse_plan, const mwSize *gauss_dims, int ny) {
	double alpha = 2.0;
	double x_sigma = 3*(5.50/6)*A;
	double y_sigma = 3*(5.50/6)*A;
	double gauss_peak_val = 3*1.61*A;
	double well_depth = 10e-3*eV;

	/// TODO: UpdateBrownianMotionGaussians
	interpolate1d_adsorbate_positions<<<1, gauss_dims[0]>>>(dev_gauss_time, dev_gauss_pos, gauss_dims[0], gauss_dims[2], dev_x0, dev_y0, t_query, ny);

	// Get the exponential of the potential
	compute_expv<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expv, expv_scale, size);

	// Apply half potential operator
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expv, size);
	cudaDeviceSynchronize();

	// Compute the forward FFT
	cufftExecZ2Z(forward_plan, dev_psi, dev_psi, CUFFT_FORWARD);

	// apply kinetic operator
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expk, size);
	cudaDeviceSynchronize();

	// Invert FFT
	cufftExecZ2Z(inverse_plan, dev_psi, dev_psi, CUFFT_INVERSE);
	complex_scale<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, 1/(double) size, size);

	/// TODO: UpdateBrownianMotionGaussians
	interpolate1d_adsorbate_positions<<<1, gauss_dims[0]>>>(dev_gauss_time, dev_gauss_pos, gauss_dims[0], gauss_dims[2], dev_x0, dev_y0, t_query, ny);

	// Get the exponential of the potential
	compute_expv<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expv, expv_scale, size);

	// Apply half potential operator
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expv, size);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Parse input parameters
	double t_query = mxGetScalar(prhs[0]);
	long long expv_ptr = mxGetScalar(prhs[1]);
	long long z_offset_ptr = mxGetScalar(prhs[2]);
	long long x0_ptr = mxGetScalar(prhs[3]);
	long long y0_ptr = mxGetScalar(prhs[4]);
	long long expk_ptr = mxGetScalar(prhs[5]);
	long long psi_ptr = mxGetScalar(prhs[6]);
	int nx = mxGetScalar(prhs[7]);
	int ny = mxGetScalar(prhs[8]);
	int nz = mxGetScalar(prhs[9]);
	int decay_type = mxGetScalar(prhs[10]);
	double A = mxGetScalar(prhs[11]);
	double eV = mxGetScalar(prhs[12]);
	double h_bar = mxGetScalar(prhs[13]);
	double dt = mxGetScalar(prhs[14]);
	double *gaussian_times = mxGetPr(prhs[15]);
	double *gaussian_positions = mxGetPr(prhs[16]);

	double expv_scale = -dt / (2 * h_bar);

	// Calculate grid size
	size_t grid_size = nx * ny * nz;

	// Get number adsorbates
	const mwSize *gauss_dims = mxGetDimensions(prhs[16]);

	// Parse the pointers
	double *dev_x0 = reinterpret_cast<double *>(x0_ptr);
	double *dev_y0 = reinterpret_cast<double *>(y0_ptr);
	double *dev_gauss_time, *dev_gauss_pos;
	cudaMallocManaged(reinterpret_cast<void **>(&dev_gauss_time), gauss_dims[0] * gauss_dims[1] * gauss_dims[2] * sizeof(double));
	cudaMallocManaged(reinterpret_cast<void **>(&dev_gauss_pos), gauss_dims[0] * gauss_dims[1] * gauss_dims[2] * sizeof(double));
	cudaMemcpy(dev_gauss_time, gaussian_times, gauss_dims[0] * gauss_dims[1] * gauss_dims[2] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_gauss_pos, gaussian_positions, gauss_dims[0] * gauss_dims[1] * gauss_dims[2] * sizeof(double), cudaMemcpyHostToDevice);

	// Plan the FFT
	cufftHandle forward_plan, inverse_plan;
	int n[3] = {nz, ny, nx};
	int idist = grid_size;
	int odist = grid_size;
	int istride = 1;
	int ostride = 1;
	int inembed[3] = {nz, ny, nx}; // MATLAB inverts rows and columns
	int onembed[3] = {nz, ny, nx};
	cufftPlanMany(&forward_plan, 3, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, 1);
	cufftPlanMany(&inverse_plan, 3, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2Z, 1);

	// TODO: PERFORM ACTUAL STEP

	cufftDestroy(forward_plan);
	cufftDestroy(inverse_plan);
	cudaFree(dev_gauss_time);
	cudaFree(dev_gauss_pos);
}

__global__ void compute_expv(myComplex *dev_expv, double scale, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		myComplex expv;
		expv.x = cos(scale * dev_expv[tid].x);
		expv.y = sin(scale * dev_expv[tid].x);

		dev_expv[tid] = expv;

		tid += blockDim.x * gridDim.x;
	}
}
