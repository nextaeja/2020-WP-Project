#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>
#include "cuda.h"

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"
#include "../Setup/cuda_setup_dynamic_potential.h"

//#include "unistd.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int NUM_GAUSSIAN_ADSORBATE_DIMENSIONS = 2;

__global__ void update_adsorbate_position(double* all_positions, double* dev_x0, double* dev_y0, int iteration, int num_adsorbates);
__device__ __host__ double _get_gaussian_adsorbate(double* data, int idx, int dim, int adsorbate, int num_adsorbates);
__global__ void compute_expv(myComplex* dev_expv, double scale, size_t size);

void split_operator_3rd_vsplit_time(myComplex* dev_psi, myComplex* dev_expv, myComplex* dev_expk,
	double* dev_x0, double* dev_y0, double* dev_z_offset, double t_query, double A, double eV,
	double expv_scale, size_t size, cufftHandle forward_plan, cufftHandle inverse_plan, const mwSize* gauss_dims,
	int nx, int ny, int nz, int decay_type, double dx, double dy, double dz, double dt,double alpha) {


	setup_dynamic_gaussian_potential(dev_expv, dev_z_offset, dev_x0, dev_y0, gauss_dims[0], nx, ny, nz, decay_type, alpha, eV, A, dx, dy, dz);

	// Get the exponential of the potential
	compute_expv << <NUM_BLOCKS, NUM_THREADS >> > (dev_expv, expv_scale, size);

	// Apply half potential operator
	complex_mul << <NUM_BLOCKS, NUM_THREADS >> > (dev_psi, dev_expv, size);

	// Compute the forward FFT
	CUDAFFT_HANDLE(cufftExecZ2Z(forward_plan, dev_psi, dev_psi, CUFFT_FORWARD));

	// apply kinetic operator
	complex_mul << <NUM_BLOCKS, NUM_THREADS >> > (dev_psi, dev_expk, size);

	// Invert FFT
	CUDAFFT_HANDLE(cufftExecZ2Z(inverse_plan, dev_psi, dev_psi, CUFFT_INVERSE));
	complex_scale << <NUM_BLOCKS, NUM_THREADS >> > (dev_psi, 1 / (double)size, size);

	/// TODO: UpdateBrownianMotionGaussians
	setup_dynamic_gaussian_potential(dev_expv, dev_z_offset, dev_x0, dev_y0, gauss_dims[0], nx, ny, nz, decay_type, alpha, eV, A, dx, dy, dz);

	// Get the exponential of the potential
	compute_expv << <NUM_BLOCKS, NUM_THREADS >> > (dev_expv, expv_scale, size);

	// Apply half potential operator
	complex_mul << <NUM_BLOCKS, NUM_THREADS >> > (dev_psi, dev_expv, size);
}

void mex_split_operator_step_3rd_vsplit_time_dependent(double t_query, myComplex* dev_expv, double* dev_z_offset, double* dev_gauss_pos, double* dev_x0, double* dev_y0, myComplex* dev_expk, myComplex* dev_psi, int nx, int ny, int nz, int decay_type, double A, double eV, double h_bar, double dt, double dx, double dy, double dz, int iteration, const mwSize* gauss_dims,double alpha) {
	// Parse input parameters
	/*double t_query = mxGetScalar(prhs[0]);
	long long expv_ptr = mxGetScalar(prhs[1]);
	long long z_offset_ptr = mxGetScalar(prhs[2]);
	long long gauss_pos_ptr = mxGetScalar(prhs[3]);
	long long x0_ptr = mxGetScalar(prhs[4]);
	long long y0_ptr = mxGetScalar(prhs[5]);
	long long expk_ptr = mxGetScalar(prhs[6]);
	long long psi_ptr = mxGetScalar(prhs[7]);
	int nx = mxGetScalar(prhs[8]);
	int ny = mxGetScalar(prhs[9]);
	int nz = mxGetScalar(prhs[10]);
	int decay_type = mxGetScalar(prhs[11]);
	double A = mxGetScalar(prhs[12]);
	double eV = mxGetScalar(prhs[13]);
	double h_bar = mxGetScalar(prhs[14]);
	double dt = mxGetScalar(prhs[15]);
	double dx = mxGetScalar(prhs[16]);
	double dy = mxGetScalar(prhs[17]);
	double dz = mxGetScalar(prhs[18]);
	int iteration = mxGetScalar(prhs[19]);
	*/
	double expv_scale = -dt / (2 * h_bar);


	// Calculate grid size
	size_t grid_size = nx * ny * nz;

	// Get number adsorbates
	//const mwSize *gauss_dims = mxGetDimensions(prhs[15]);

	// Parse the pointers
	/*myComplex *dev_expv = reinterpret_cast<myComplex *>(expv_ptr);
	double *dev_z_offset = reinterpret_cast<double *>(z_offset_ptr);
	double *dev_gauss_pos = reinterpret_cast<double *>(gauss_pos_ptr);
	double *dev_x0 = reinterpret_cast<double *>(x0_ptr);
	double *dev_y0 = reinterpret_cast<double *>(y0_ptr);
	myComplex *dev_expk = reinterpret_cast<myComplex *>(expk_ptr);
	myComplex *dev_psi = reinterpret_cast<myComplex *>(psi_ptr);
	*/
	// Plan the FFT
	cufftHandle forward_plan, inverse_plan;
	int n[3] = { nz, ny, nx };
	int idist = grid_size;
	int odist = grid_size;
	int istride = 1;
	int ostride = 1;
	int inembed[3] = { nz, ny, nx }; // MATLAB inverts rows and columns
	int onembed[3] = { nz, ny, nx };
	CUDAFFT_HANDLE(cufftPlanMany(&forward_plan, 3, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, 1));
	CUDAFFT_HANDLE(cufftPlanMany(&inverse_plan, 3, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2Z, 1));

	// Compute the x and y positions of the adsorbates
	update_adsorbate_position << <1, gauss_dims[0] >> > (dev_gauss_pos, dev_x0, dev_y0, iteration, gauss_dims[0]);

	split_operator_3rd_vsplit_time(dev_psi, dev_expv, dev_expk, dev_x0, dev_y0, dev_z_offset, t_query, A, eV, expv_scale, grid_size, forward_plan, inverse_plan, gauss_dims, nx, ny, nz, decay_type, dx, dy, dz, dt,alpha);

	CUDAFFT_HANDLE(cufftDestroy(forward_plan));
	CUDAFFT_HANDLE(cufftDestroy(inverse_plan));
}

__global__ void update_adsorbate_position(double* all_positions, double* dev_x0, double* dev_y0, int iteration, int num_adsorbates) {
	int adsorbate = blockIdx.x * blockDim.x + threadIdx.x;

	while (adsorbate < num_adsorbates) {
		dev_x0[adsorbate] = _get_gaussian_adsorbate(all_positions, iteration, 0, adsorbate, num_adsorbates);
		dev_y0[adsorbate] = _get_gaussian_adsorbate(all_positions, iteration, 1, adsorbate, num_adsorbates);

		adsorbate += blockDim.x * gridDim.x;
	}
}

// The gaussian position array is a 3D one, return the correct value for a given
//	adsorbate number and dimension (x or y)
__device__ __host__ double _get_gaussian_adsorbate(double* data, int idx, int dim, int adsorbate, int num_adsorbates) {
	if (idx < 0) {
		return -1.0;
	}

	int tot_idx = adsorbate + dim * num_adsorbates + idx * num_adsorbates * NUM_GAUSSIAN_ADSORBATE_DIMENSIONS;

	return data[tot_idx];
}

__global__ void compute_expv(myComplex* dev_expv, double scale, size_t size) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	while (tid < size) {
		myComplex expv;
		expv.x = cos(scale * dev_expv[tid].x);
		expv.y = sin(scale * dev_expv[tid].x);

		dev_expv[tid] = expv;

		tid += blockDim.x * gridDim.x;
	}
}






void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	char cwdtmp[256];

	//getcwd(cwd, sizeof(cwd));
	char* cwd = "C:\\Users\\jackl\\Documents\\MATLAB\\summer 2020\\2020-WP-Project-6_optimizing_setup_dynamic_gaussian_potential\\Beta4_2\\mexcudawhile";//make general
	char* l = strrchr(cwd, '\\');
	cwd[strlen(cwd) - strlen(l) + 1] = '\0';

	strcpy(cwdtmp, cwd);

	FILE* setup = fopen(strcat(cwdtmp, "setup.txt"), "r");//navigates up to setup file


	char buff[1024];
	
	ps = e-12;

	fgets(buff, sizeof(buff), setup);
	expv_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	z_offset_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	gauss_pos_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	x0_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	y0_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	expk_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	psi_ptr = strtoll(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	nx = atoi(buff);

	fgets(buff, sizeof(buff), setup);
	ny = atoi(buff);

	fgets(buff, sizeof(buff), setup);
	nz = atoi(buff);
	
	fgets(buff, sizeof(buff), setup);
	decay_type = atoi(buff);

	fgets(buff, sizeof(buff), setup);
	A = strtod(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	eV = strtod(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	h_bar = strtod(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	dt = strtod(buff, NULL) * ps;

	fgets(buff, sizeof(buff), setup);
	dx = strtod(buff, NULL) * A;

	fgets(buff, sizeof(buff), setup);
	dy = strtod(buff, NULL) * A;

	fgets(buff, sizeof(buff), setup);
	dz = strtod(buff, NULL) * A;

	fgets(buff, sizeof(buff), setup);
	gfxSteps = strtol(buff, NULL);

	fgets(buff, sizeof(buff), setup);
	t = strtod(buff, NULL) * ps;

	fgets(buff, sizeof(buff), setup);
	alpha = strtod(buff, NULL);

	fclose(setup);



	int it = 0;
	//int* sum;

	myComplex* dev_expv = reinterpret_cast<myComplex*>(expv_ptr);
	double* dev_z_offset = reinterpret_cast<double*>(z_offset_ptr);
	double* dev_gauss_pos = reinterpret_cast<double*>(gauss_pos_ptr);
	double* dev_x0 = reinterpret_cast<double*>(x0_ptr);
	double* dev_y0 = reinterpret_cast<double*>(y0_ptr);
	myComplex* dev_expk = reinterpret_cast<myComplex*>(expk_ptr);
	myComplex* dev_psi = reinterpret_cast<myComplex*>(psi_ptr);

	//const mwSize* gauss_dims = mxGetDimensions(prhs[15]);

	const size_t gauss_dims = sizeof(dx);//no idea if this is actually a reasonable replacement

	while (it <= gfxSteps) {
		/*double totprob = 0;
		for (int x = 0;x < nx;x++) {
			for (int y = 0;y < ny;y++) {
				for (int z = 0;z < nz;z++) {
					for (int im = 0;im < 2;im++) {
						totprob = totprob + psi[x][y][z][im] * psi[x][y][z][im];
					}
				}
			}
		}


		if (!(fabs(totprob - 1) < eps)) {
			nrhs = 1;
			prhs[0] = 5;//make actual error
		}*/

		//* sum = 0;



		mex_split_operator_step_3rd_vsplit_time_dependent(t, dev_expv, dev_z_offset, dev_gauss_pos, dev_x0, dev_y0, dev_expk, dev_psi, nx, ny, nz, decay_type, A, eV, h_bar, dt, dx, dy, dz, it, gauss_dims,alpha);
		it++;
		t = t + dt;
	}

}