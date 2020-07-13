#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>     // Needed for memcpy()


// Evaluate the morse potential for a single value of z
__device__ inline double morse_potential(double z, double wellDepth, double wellMinZPt, double a) {
    return wellDepth * (pow(1 - exp(-a*(z-wellMinZPt)), 2) - 1);
}

__device__ inline double gaussian_fun(double x, double y, double ads_x, double ads_y, double x_sigma, double y_sigma) {
    double x_arg = - pow(x - ads_x, 2) / (2 * pow(x_sigma, 2));
    double y_arg = - pow(y - ads_y, 2) / (2 * pow(y_sigma, 2));

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
        potential[tid] = 0.0;
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
    double *x0 = mxGetPr(prhs[6]), *dev_x0;
    double *y0 = mxGetPr(prhs[7]), *dev_y0;
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

    // Extract the number of adsorbates
    const int adsorbate_num = mxGetM(prhs[6]);

    // Copy the adsorbate position in GPU
    cudaMalloc((void**) &dev_x0, adsorbate_num * sizeof(double));
    cudaMalloc((void**) &dev_y0, adsorbate_num * sizeof(double));
    cudaMemcpy((void *) dev_x0, (void *) x0, adsorbate_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void *) dev_y0, (void *) y0, adsorbate_num * sizeof(double), cudaMemcpyHostToDevice);

    // Setup the number of dimensions in the output array
    const mwSize potential_dims[] = {nx, ny, nz};
    size_t data_length = nx * ny * nz;

    // Allocate a block of memory in the GPU to save the potential in
    double *potential;
    cudaMalloc((void**) &potential, data_length * sizeof(double));
    // TODO add an error check for bad allocation
    // TODO check out cudaMalloc3D

    // Allocate a 2D grid for z offsets and xy grid
    double *z_offset;
    cudaMalloc((void**) &z_offset, nx * ny * sizeof(double));

    // Initialize potential to zero
    //initialize_potential<<<128, 128>>>(potential, data_length);

    // Compute the z offset
    dim3 dimBlock(1, 1);
    dim3 dimGrid_xy(nx, ny);
    compute_z_offset<<<dimGrid_xy, dimBlock>>>(z_offset, nx, ny, dx, dy, dev_x0, dev_y0, adsorbate_num, xSigma, ySigma, gaussPeakVal);

    // Compute the 3d potential
    dim3 dimGrid_3d(nx*ny, nz);
    compute_potential<<<dimGrid_3d, dimBlock>>>(potential, z_offset, nx, ny, nz, dz, wellDepth, wellMinZPt, a);

    // Create a 3D array of doubles of dimensions [nx, ny, nz] for the computed potential
    double *final_potential;
    plhs[0] = mxCreateNumericArray(3, potential_dims, mxDOUBLE_CLASS, mxREAL);
    final_potential = (double *) mxGetPr(plhs[0]);

    // Copy the computed array into the MATLAB one
    cudaMemcpy(final_potential, potential, data_length * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(potential);
    cudaFree(z_offset);
    cudaFree((void *) dev_x0);
    cudaFree((void *) dev_y0);
}
