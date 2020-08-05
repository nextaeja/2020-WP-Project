#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"
#include "../MEX_helpers/interpolation1d.h"
#include "../Setup/cuda_setup_dynamic_potential.h"

#define NDIMS 3

__global__ void update_adsorbate_position(double *all_positions, double *dev_x0, double *dev_y0, int iteration, int num_adsorbates);
__device__ __host__ double _get_gaussian_adsorbate(double *data, int idx, int dim, int adsorbate, int num_adsorbates);
__global__ void compute_expv(myComplex *dev_expv, double scale, size_t size);

void split_operator_3rd_vsplit_time(myComplex *dev_psi, myComplex *dev_expv, myComplex *dev_expk,
		double *dev_x0, double *dev_y0, double *dev_z_offset, double t_query, double A, double eV,
		double expv_scale, size_t size, cufftHandle forward_plan, cufftHandle inverse_plan, const mwSize *gauss_dims,
		int nx, int ny, int nz, int decay_type, double dx, double dy, double dz, double dt) {
	double alpha = 2.0;

	/// TODO: UpdateBrownianMotionGaussians
	setup_dynamic_gaussian_potential(dev_expv, dev_z_offset, dev_x0, dev_y0, gauss_dims[0], nx, ny, nz, decay_type, alpha, eV, A, dx, dy, dz);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// Get the exponential of the potential
	compute_expv<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expv, expv_scale, size);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// Apply half potential operator
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expv, size);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// Compute the forward FFT
	CUDAFFT_HANDLE(cufftExecZ2Z(forward_plan, dev_psi, dev_psi, CUFFT_FORWARD));
	CUDA_HANDLE(cudaDeviceSynchronize());

	// apply kinetic operator
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expk, size);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// Invert FFT
	CUDAFFT_HANDLE(cufftExecZ2Z(inverse_plan, dev_psi, dev_psi, CUFFT_INVERSE));
	CUDA_HANDLE(cudaDeviceSynchronize());
	complex_scale<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, 1/(double) size, size);
	CUDA_HANDLE(cudaDeviceSynchronize());

	/// TODO: UpdateBrownianMotionGaussians
	setup_dynamic_gaussian_potential(dev_expv, dev_z_offset, dev_x0, dev_y0, gauss_dims[0], nx, ny, nz, decay_type, alpha, eV, A, dx, dy, dz);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// Get the exponential of the potential
	compute_expv<<<NUM_BLOCKS, NUM_THREADS>>>(dev_expv, expv_scale, size);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// Apply half potential operator
	complex_mul<<<NUM_BLOCKS, NUM_THREADS>>>(dev_psi, dev_expv, size);
	CUDA_HANDLE(cudaDeviceSynchronize());
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
	double *gaussian_positions = mxGetPr(prhs[15]);
	double dx = mxGetScalar(prhs[16]);
	double dy = mxGetScalar(prhs[17]);
	double dz = mxGetScalar(prhs[18]);
	int iteration = mxGetScalar(prhs[19]);

	double expv_scale = -dt / (2 * h_bar);

	size_t total, free;
	CUDA_HANDLE(cudaMemGetInfo(&free, &total));
	//mexPrintf("Free: %.2lf of %.2lf (%.2lf%%)\n", free / (double) (1 << 30), total / (double) (1 << 30), free / (double) total * 100.0);

	// Calculate grid size
	size_t grid_size = nx * ny * nz;

	// Get number adsorbates
	const mwSize *gauss_dims = mxGetDimensions(prhs[15]);

	// Parse the pointers
	myComplex *dev_expv = reinterpret_cast<myComplex *>(expv_ptr);
	double *dev_z_offset = reinterpret_cast<double *>(z_offset_ptr);
	double *dev_x0 = reinterpret_cast<double *>(x0_ptr);
	double *dev_y0 = reinterpret_cast<double *>(y0_ptr);
	myComplex *dev_expk = reinterpret_cast<myComplex *>(expk_ptr);
	myComplex *dev_psi = reinterpret_cast<myComplex *>(psi_ptr);

	// Copy the adsorbate position and times into GPU
	// TODO: move this out of this function. Only execute at beginning
	double *dev_gauss_pos;
	CUDA_HANDLE(cudaMallocManaged(reinterpret_cast<void **>(&dev_gauss_pos), gauss_dims[0] * gauss_dims[1] * gauss_dims[2] * sizeof(double)));
	CUDA_HANDLE(cudaMemcpy(dev_gauss_pos, gaussian_positions, gauss_dims[0] * gauss_dims[1] * gauss_dims[2] * sizeof(double), cudaMemcpyHostToDevice));

	// Plan the FFT
	cufftHandle forward_plan, inverse_plan;
	int n[3] = {nz, ny, nx};
	int idist = grid_size;
	int odist = grid_size;
	int istride = 1;
	int ostride = 1;
	int inembed[3] = {nz, ny, nx}; // MATLAB inverts rows and columns
	int onembed[3] = {nz, ny, nx};
	CUDAFFT_HANDLE(cufftPlanMany(&forward_plan, 3, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, 1));
	CUDAFFT_HANDLE(cufftPlanMany(&inverse_plan, 3, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2Z, 1));

	// Compute the x and y positions of the adsorbates
	update_adsorbate_position<<<1, gauss_dims[0]>>>(dev_gauss_pos, dev_x0, dev_y0, iteration, gauss_dims[0]);
	CUDA_HANDLE(cudaDeviceSynchronize());

	// TODO: PERFORM ACTUAL STEP
	split_operator_3rd_vsplit_time(dev_psi, dev_expv, dev_expk, dev_x0, dev_y0, dev_z_offset, t_query, A, eV, expv_scale, grid_size, forward_plan, inverse_plan, gauss_dims, nx, ny, nz, decay_type, dx, dy, dz, dt);

	CUDAFFT_HANDLE(cufftDestroy(forward_plan));
	CUDAFFT_HANDLE(cufftDestroy(inverse_plan));
	CUDA_HANDLE(cudaFree(dev_gauss_pos));
}

__global__ void update_adsorbate_position(double *all_positions, double *dev_x0, double *dev_y0, int iteration, int num_adsorbates) {
	int adsorbate = blockIdx.x * blockDim.x + threadIdx.x;

	while (adsorbate < num_adsorbates) {
		dev_x0[adsorbate] = _get_gaussian_adsorbate(all_positions, iteration, 0, adsorbate, num_adsorbates);
		dev_y0[adsorbate] = _get_gaussian_adsorbate(all_positions, iteration, 1, adsorbate, num_adsorbates);

		adsorbate += blockDim.x * gridDim.x;
	}
}

// The gaussian position array is a 3D one, return the correct value for a given
//	adsorbate number and dimension (x or y)
__device__ __host__ double _get_gaussian_adsorbate(double *data, int idx, int dim, int adsorbate, int num_adsorbates) {
	if (idx < 0) {
		return -1.0;
	}

	int tot_idx = adsorbate + dim*num_adsorbates + idx*num_adsorbates*NUM_GAUSSIAN_ADSORBATE_DIMENSIONS;

	return data[tot_idx];
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
