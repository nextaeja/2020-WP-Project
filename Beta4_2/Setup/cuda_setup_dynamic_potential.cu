/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "../MEX_helpers/cuda_helper.h"
#include "../MEX_helpers/complex.h"
#include "cuda_setup_dynamic_potential.h"

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

    myComplex point_potential;
    point_potential.x = morse_potential(z, wellDepth, wellMinZPt, a);
    point_potential.y = 0.0;

    potential[idx] = point_potential;
}

/*
 *  Entry point for the C version of SetupDynamicGaussianPotential
 */
void setup_dynamic_gaussian_potential(myComplex *dev_potential, double *dev_z_offset, double *dev_x0, double *dev_y0,
        int adsorbate_num, int nx, int ny, int nz, int decayType, int inParameterIfNeeded, double eV, double A,
        double dx, double dy, double dz) {

    // Setup helper variables
    double wellMinZPt = 2*A;
    double zCharacteristic = (1/2.06)*A;
    double wellDepth = 10e-3*eV;
    double gaussPeakVal = 3*1.61*A;
    double xSigma = pow(3*(5.50/6)*A, 2.0);
    double ySigma = pow(3*(5.50/6)*A, 2.0);
    double Vmax = wellDepth*exp(wellMinZPt/zCharacteristic);
    double a = (1/wellMinZPt)*log(1+sqrt(1+(Vmax/wellDepth)));

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
