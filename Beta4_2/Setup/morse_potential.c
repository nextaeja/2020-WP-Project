#include <math.h>
#include <matrix.h>
#include <mex.h>


/*
 *  Function to evaluate the morse potential on a 3D grid.
 *  INPUTS:
 *   - zEffective -> 3D array of effective z values
 *   - wellDepth  -> float
 *   - wellMinZPt -> float
 *   - a          -> float
 *
 *  OUTPUTS:
 *   - zEffective3D
 */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    // TODO: add parameter number check

    mxArray *z_array;
    double *z_values, well_depth, well_min_zpt, a;
    const mwSize *grid_size;
    
    // Get the potential parameters
    well_depth = mxGetScalar(prhs[1]);
    well_min_zpt = mxGetScalar(prhs[2]);
    a = mxGetScalar(prhs[3]);
    
    //mexPrintf("%d %d\n", nrhs, nlhs);
    //mexPrintf("%e %e %e\n\n", well_depth, well_min_zpt, a);
    
    grid_size = mxGetDimensions(prhs[0]);
    mexPrintf("%d %d %d\n", (int) grid_size[0], (int) grid_size[1], (int) grid_size[2]);

    return;
    
    // Get the array of z
    z_array = mxDuplicateArray(prhs[0]);
    z_values = mxGetPr(z_array);
}




