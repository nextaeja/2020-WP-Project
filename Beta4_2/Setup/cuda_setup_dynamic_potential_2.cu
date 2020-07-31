#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "../MEX_helpers/cuda_helper.h"
#include "../MEX_helpers/complex.h"

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
__device__ inline double compute_single_offset(double x, double y, double gaussian_peak_value,
        double *x0, double *y0, int adsorbate_num, double x_sigma, double y_sigma) {
    double z_offset = 0.0;

    // Loop over each adsorbate
    for (int i=0; i<adsorbate_num; i++) {
        z_offset += gaussian_peak_value * gaussian_fun(x, y, x0[i], y0[i], x_sigma, y_sigma);
    }

    return z_offset;
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

__global__ void compute_potential(myComplex *potential, double *z_offset, size_t nxy, size_t nz,
        double dz, double wellDepth, double wellMinZPt, double a) {
    // Calculate the coordinate of the specific thread
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row * nxy + col;
    if (row >= nz || col >= nxy || idx >= nxy*nz)
        return;

    double z = dz * row - z_offset[col];

    double point_potential = morse_potential(z, wellDepth, wellMinZPt, a);

    myComplex w;
    w.x = point_potential;
    w.y = 0.0;
    potential[idx] = w;
}

/*
 *  Entry point for the C version of SetupDynamicGaussianPotential
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Parse the parameters
    long long potential_ptr = mxGetScalar(prhs[0]);
    long long z_offset_ptr = mxGetScalar(prhs[1]);
    long long x0_ptr = mxGetScalar(prhs[2]);
    long long y0_ptr = mxGetScalar(prhs[3]);
    int decayType = mxGetScalar(prhs[4]);
    int inParameterIfNeeded = mxGetScalar(prhs[5]);
    double xSigma = mxGetScalar(prhs[6]);
    double ySigma = mxGetScalar(prhs[7]);
    double gaussPeakVal = mxGetScalar(prhs[8]);
    double wellDepth = mxGetScalar(prhs[9]);
    double *x0 = mxGetPr(prhs[10]);
    double *y0 = mxGetPr(prhs[11]);
    double dx = mxGetScalar(prhs[12]);
    double dy = mxGetScalar(prhs[13]);
    double dz = mxGetScalar(prhs[14]);
    double A = mxGetScalar(prhs[15]);
    size_t nx = mxGetScalar(prhs[16]);
    size_t ny = mxGetScalar(prhs[17]);
    size_t nz = mxGetScalar(prhs[18]);

    // Get the pointer to the GPU arrays allocated earlier
    myComplex *dev_potential = reinterpret_cast<myComplex *>(potential_ptr);
    double *dev_z_offset = reinterpret_cast<double *>(z_offset_ptr);
    double *dev_x0 = reinterpret_cast<double *>(x0_ptr);
    double *dev_y0 = reinterpret_cast<double *>(y0_ptr);

    // Setup helper variables
    double wellMinZPt = 2*A;
    double zCharacteristic = (1/2.06)*A;
    double Vmax = wellDepth*exp(wellMinZPt/zCharacteristic);
    double a = (1/wellMinZPt)*log(1+sqrt(1+(Vmax/wellDepth)));
    xSigma = pow(xSigma, 2.0);
    ySigma = pow(ySigma, 2.0);

    // Extract the number of adsorbates
    const int adsorbate_num = mxGetM(prhs[10]);

    // Copy the adsorbate positions in GPU
    cudaMemcpy(dev_x0, x0, adsorbate_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y0, y0, adsorbate_num * sizeof(double), cudaMemcpyHostToDevice);

    // Compute the z offset
    int blocks_x = (nx + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    int blocks_y = (ny + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    dim3 dimBlock_xy(blocks_x, blocks_y);
    dim3 dimThread(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    compute_z_offset<<<dimBlock_xy, dimThread>>>(dev_z_offset, nx, ny, dx, dy, dev_x0, dev_y0, adsorbate_num, xSigma, ySigma, gaussPeakVal);

    // Compute the 3d potential
    int blocks_xy = (nx*ny + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    int blocks_z = (nz + THREADS_PER_BLOCK - 1) /  THREADS_PER_BLOCK;
    dim3 dimBlock_3d(blocks_xy, blocks_z);
    compute_potential<<<dimBlock_3d, dimThread>>>(dev_potential, dev_z_offset, nx * ny, nz, dz, wellDepth, wellMinZPt, a);
}
