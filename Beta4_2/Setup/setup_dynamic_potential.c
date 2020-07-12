#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>     // Needed for memcpy()

const int REQUIRED_INPUT_ARGUMENTS = 15;

// Perform a sanity check on the input data
void error_check(int nlhs, int nrhs, const mxArray *prhs[]) {
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:lhs",
                          "The function needs exactly one output variable");
    }

    if (nrhs < REQUIRED_INPUT_ARGUMENTS) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:rhs",
                          "Too few arguments given");
    }

    if (nrhs > REQUIRED_INPUT_ARGUMENTS) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:rhs",
                          "Too many arguments given");
    }

    // Check that x0 is a matrix of size [M 1]
    if (mxGetN(prhs[6]) != 1 && mxGetM(prhs[6]) > 0) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:MalformedInput",
                          "Adsorbate x position needs to be a column vector");
    }

    // Check that y0 is a matrix of size [M 1]
    if (mxGetN(prhs[7]) != 1 && mxGetM(prhs[7]) > 0) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:MalformedInput",
                          "Adsorbate y position needs to be a column vector");
    }

    // Check that x0 and y0 have the same dimension
    if (mxGetM(prhs[6]) != mxGetM(prhs[7])) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:MalformedInput",
                          "Adsorbate x and y position arrays need to be of the same size");
    }
}

// Evaluate the morse potential for a single value of z
double morse_potential(double z, const double wellDepth, const double wellMinZPt, const double a) {
    return wellDepth * (pow(1 - exp(-a*(z-wellMinZPt)), 2) - 1);
}

double gaussian_fun(double x, double y, double ads_x, double ads_y, double x_sigma, double y_sigma) {
    double x_arg = - pow(x - ads_x, 2) / (2 * pow(x_sigma, 2));
    double y_arg = - pow(y - ads_y, 2) / (2 * pow(y_sigma, 2));

    return exp(x_arg + y_arg);
}

// Compute the z offest from adsorbates
double compute_z_offset(double x, double y, const double gaussian_peak_value,
        const double *adsorbate_x, const double *adsorbate_y, const int adsorbate_num,
        const double x_sigma, const double y_sigma) {
    double z_offset = 0.0;

    // Loop over each adsorbate
    for (int i=0; i<adsorbate_num; i++) {
        z_offset += gaussian_peak_value * gaussian_fun(x, y, adsorbate_x[i], adsorbate_y[i], x_sigma, y_sigma);
    }

    return z_offset;
}

/*
 *  Entry point for the C version of SetupDynamicGaussianPotential
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double wellMinZPt, zCharacteristic, Vmax, a;

    // Check the correctness of the number of inputs and outputs
    error_check(nlhs, nrhs, prhs);

    // Parse the parameters
    const int decayType = mxGetScalar(prhs[0]);
    const int inParameterIfNeeded = mxGetScalar(prhs[1]);
    const double xSigma = mxGetScalar(prhs[2]);
    const double ySigma = mxGetScalar(prhs[3]);
    const double gaussPeakVal = mxGetScalar(prhs[4]);
    const double wellDepth = mxGetScalar(prhs[5]);
    const double *x0 = mxGetPr(prhs[6]);
    const double *y0 = mxGetPr(prhs[7]);
    const double dx = mxGetScalar(prhs[8]);
    const double dy = mxGetScalar(prhs[9]);
    const double dz = mxGetScalar(prhs[10]);
    const double A = mxGetScalar(prhs[11]);
    const int nx = mxGetScalar(prhs[12]);
    const int ny = mxGetScalar(prhs[13]);
    const int nz = mxGetScalar(prhs[14]);

    // Extract the number of adsorbates
    const int adsorbate_num = mxGetM(prhs[6]);

    // Setup the number of dimensions in the output array
    const mwSize potential_dims[] = {nx, ny, nz};
    int data_length = nx * ny * nz;

    // Allocate a block of memory to save the potential in
    double *potential = (double *) malloc(data_length * sizeof(double));
    if (potential == NULL) {
        mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:MemoryFailure",
                          "Failure to allocate memory");
    }

    // Initialize potential to zero
    for (int i=0; i<data_length; i++) {
        potential[i] = 0;
    }

    // Change behaviour depending of decay type
    switch (decayType) {
        case 1:
            mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:NotImplemented",
                              "This decay type has not been implemented yet");
            break;

        // Morse potential
        case 2:
            // Setup helper variables
            wellMinZPt = 2*A;
            zCharacteristic = (1/2.06)*A;
            Vmax = wellDepth*exp(wellMinZPt/zCharacteristic);
            a = (1/wellMinZPt)*log(1+sqrt(1+(Vmax/wellDepth)));

            // For each space position compute the offset
            for (int i=0; i<nx; i++) {
                double x = i * dx;

                for (int j=0; j<ny; j++) {
                    double y = j * dy;
                    double offset = compute_z_offset(x, y, gaussPeakVal, x0, y0, adsorbate_num, xSigma, ySigma);

                    for (int k=0; k<nz; k++) {
                        double z = k * dz;
                        double eff_z = z;

                        // Go from 3D index to linear index
                        int lin_idx = i + j*nx + k*nx*ny;

                        // Compute the effective value of z from the adsorbate positions
                        eff_z -= offset;

                        // Calculate the potential at the point
                        potential[lin_idx] = morse_potential(eff_z, wellDepth, wellMinZPt, a);
                    }
                }
            }

            break;

        case 3:
            mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:NotImplemented",
                              "This decay type has not been implemented yet");
            break;

        // Invalid decay type
        default:
            mexErrMsgIdAndTxt("MATLAB:SetupDynamicGaussianPotential:InvalidDecay",
                              "Invalid decay type given");
    }

    // Create a 3D array of doubles of dimensions [nx, ny, nz] for the computed potential
    double *final_potential;
    plhs[0] = mxCreateNumericArray(3, potential_dims, mxDOUBLE_CLASS, mxREAL);
    final_potential = (double *) mxGetPr(plhs[0]);

    // Copy the computed array into the MATLAB one
    memcpy(final_potential, potential, data_length * sizeof(double));
    free(potential);
}
