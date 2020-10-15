/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#include <mex.h>

#include "helper.h"
#include "../MEX_helpers/cuda_helper.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get how many C and CUDA arrays are to be freed
    int num_standard = mxGetScalar(prhs[0]);
    int num_cuda = mxGetScalar(prhs[1]);

    double *standard_arrays = mxGetPr(prhs[2]);
    double *cuda_arrays = mxGetPr(prhs[3]);

    // Free all the standard arrays
    for (int i=0; i<num_standard; i++) {
        // Convert the int value into a pointer
        double *pointer = int_to_pointer(standard_arrays[i]);

        // Free the pointer
        free(pointer);
    }

    // Free all the cuda arrays
    for (int i=0; i<num_cuda; i++) {
        // Convert the int value into a pointer
        double *pointer = int_to_pointer(cuda_arrays[i]);

        // Free the pointer
        CUDA_HANDLE(cudaFree(pointer));
    }
}
