#define LINUX_PROFILE

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>     // Needed for memcpy()

#include "helper.h"
#include "../MEX_helpers/cuda_helper.h"

#ifdef LINUX_PROFILE
#include <time.h>  // Only for linux, profiling purposes
#endif

// Evaluate the morse potential for a single value of z
__device__ inline double morse_potential(double z, double wellDepth, double wellMinZPt, double a) {
    return wellDepth * (pow(1 - exp(-a*(z-wellMinZPt)), 2) - 1);
}

__device__ inline double gaussian_fun(double x, double y, double ads_x, double ads_y, double x_sigma, double y_sigma) {
    double x_arg = - pow(x - ads_x, 2) / (2 * x_sigma);
    double y_arg = - pow(y - ads_y, 2) / (2 * y_sigma);

    return exp(x_arg + y_arg);
}

// Compute the z offest from adsorbates
__device__ double compute_single_offset(double x, double y, double gaussian_peak_value,
        double *x0, double *y0, int adsorbate_num, double x_sigma, double y_sigma) {
    double z_offset = 0.0;

    // Loop over each adsorbate
    for (int i=0; i<adsorbate_num; i++) {
        z_offset += gaussian_peak_value * gaussian_fun(x, y, x0[i], y0[i], x_sigma, y_sigma);
    }

    return z_offset;
}

__global__ void initialize_potential(double *potential, size_t data_length) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < data_length) {
        potential[tid] = 42.0;
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void compute_z_offset(double *offsets, size_t nx, size_t ny, double dx, double dy,
        double *x0, double *y0, int adsorbate_num, double xSigma, double ySigma,
        double gaussPeakVal) {
    // Calculate the coordinate of the specific thread
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int xy_idx = row * nx + col;
    if (row >= ny || col >= nx)
        return;

    // Get coordinates
    double x = col * dx;
    double y = row * dy;

    double offset = compute_single_offset(x, y, gaussPeakVal, x0, y0, adsorbate_num, xSigma, ySigma);

    offsets[xy_idx] = offset;
}

__global__ void compute_potential(double *potential, double *z_offset, size_t nx, size_t ny, size_t nz,
        double dz, double wellDepth, double wellMinZPt, double a) {
    // Calculate the coordinate of the specific thread
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row * nx * ny + col;
    if (row >= nz || col >= nx*ny || idx >= nx*ny*nz)
        return;

    double z = dz * row - z_offset[col];

    potential[idx] = morse_potential(z, wellDepth, wellMinZPt, a);
}

#ifdef LINUX_PROFILE
struct timespec start, end;
void start_time() {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
}

double us_elapsed() {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double time_taken;
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    return (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-6;
}
#endif

/*
 *  Entry point for the C version of SetupDynamicGaussianPotential
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Parse the parameters
    int decayType = mxGetScalar(prhs[0]);
    int inParameterIfNeeded = mxGetScalar(prhs[1]);
    double xSigma = mxGetScalar(prhs[2]);
    double ySigma = mxGetScalar(prhs[3]);
    double gaussPeakVal = mxGetScalar(prhs[4]);
    double wellDepth = mxGetScalar(prhs[5]);
    double *x0 = mxGetPr(prhs[6]);
    double *y0 = mxGetPr(prhs[7]);
    double dx = mxGetScalar(prhs[8]);
    double dy = mxGetScalar(prhs[9]);
    double dz = mxGetScalar(prhs[10]);
    double A = mxGetScalar(prhs[11]);
    size_t nx = mxGetScalar(prhs[12]);
    size_t ny = mxGetScalar(prhs[13]);
    size_t nz = mxGetScalar(prhs[14]);

    // Setup helper variables
    double wellMinZPt = 2*A;
    double zCharacteristic = (1/2.06)*A;
    double Vmax = wellDepth*exp(wellMinZPt/zCharacteristic);
    double a = (1/wellMinZPt)*log(1+sqrt(1+(Vmax/wellDepth)));
    xSigma = pow(xSigma, 2.0);
    ySigma = pow(ySigma, 2.0);

    // Extract the number of adsorbates
    const int adsorbate_num = mxGetM(prhs[6]);

    // Setup the number of dimensions in the output array
    const mwSize potential_dims[] = {nx, ny, nz};
    size_t data_length = nx * ny * nz;

    /*
     * Pointers to the 3d array to store the potential in and the z offsets
     * Declared static so the memory is only allocated once every MATLAB session
     *  instead of every function execution. This should reduce the runtime.
     * The algorithm overwrites every cell in the potential, so initializing to zero
     *  at startup is not needed.
     */
    static double *potential=NULL;
    static double *z_offset=NULL;
    static double *dev_x0=NULL;
    static double *dev_y0=NULL;

    // If running for the first time in session allocate the array
    if (potential == NULL) {
        cudaMallocManaged(&potential, data_length * sizeof(double));
        // TODO add an error check for bad allocation
    }
    if (z_offset == NULL) {
        cudaMallocManaged((void**) &z_offset, nx * ny * sizeof(double));
    }
    if (dev_x0 == NULL) {
        cudaMallocManaged((void**) &dev_x0, adsorbate_num * sizeof(double));
    }
    if (dev_y0 == NULL) {
        cudaMallocManaged((void**) &dev_y0, adsorbate_num * sizeof(double));
    }

    // Copy the adsorbate position in GPU
    cudaMemcpy((void *) dev_x0, (void *) x0, adsorbate_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void *) dev_y0, (void *) y0, adsorbate_num * sizeof(double), cudaMemcpyHostToDevice);


#ifdef LINUX_PROFILE
    double *timing;
    const mwSize time_size[] = {6};
    plhs[1] = mxCreateNumericArray(1, time_size, mxDOUBLE_CLASS, mxREAL);
    timing = (double *) mxGetPr(plhs[1]);
#endif

#ifdef LINUX_PROFILE
    start_time();
#endif
#ifdef LINUX_PROFILE
    timing[0] = us_elapsed();
#endif

    // Compute the z offset
    int blocks_x = (nx + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    int blocks_y = (ny + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    dim3 dimBlock_xy(blocks_x, blocks_y);
    dim3 dimThread(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    compute_z_offset<<<dimBlock_xy, dimThread>>>(z_offset, nx, ny, dx, dy, dev_x0, dev_y0, adsorbate_num, xSigma, ySigma, gaussPeakVal);

    // Compute the 3d potential
    int blocks_xy = (nx*ny + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    int blocks_z = (nz + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    dim3 dimBlock_3d(blocks_xy, blocks_z);
    compute_potential<<<dimBlock_3d, dimThread>>>(potential, z_offset, nx, ny, nz, dz, wellDepth, wellMinZPt, a);

    // Create a 3D array of doubles of dimensions [nx, ny, nz] for the computed potential
    double *final_potential;
    plhs[0] = mxCreateNumericArray(3, potential_dims, mxDOUBLE_CLASS, mxREAL);
    final_potential = (double *) mxGetPr(plhs[0]);

    // Copy the computed array into the MATLAB one
    cudaMemcpy(final_potential, potential, data_length * sizeof(double), cudaMemcpyDeviceToHost);

    // Create output variables to store the pointers to the potentials and avoid memory leaks
    plhs[1] = mxCreateNumericMatrix(1, 4, mxINT64_CLASS, mxREAL);
    long long *all_potentials = (long long *) mxGetPr(plhs[1]);

    // Set the union to the value of the pointer and the return array to the correesponding
    // integer value.
    all_potentials[0] = pointer_to_int(potential);
    all_potentials[1] = pointer_to_int(z_offset);
    all_potentials[2] = pointer_to_int(dev_x0);
    all_potentials[3] = pointer_to_int(dev_y0);
}
