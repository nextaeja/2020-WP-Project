#ifndef __INTERPOLATION_1D__
#define __INTERPOLATION_1D__

#include <mex.h>

const int NUM_GAUSSIAN_ADSORBATE_DIMENSIONS = 2;

__global__ void interpolate1d_adsorbate_positions(double *gaussian_times, double *gaussian_positions, int num_adsorbates, int num_gaussian_positions, double *dev_x0, double *dev_y0, double t_query, int ny);
__device__ __host__ double get_gaussian_adsorbate(double *data, int idx, int dim, int adsorbate, int num_adsorbates);

__device__ __host__ int lin_left_locate_3d(double *data, int size, double to_locate, int num_adsorbates, int adsorbate, int dim);
__device__ __host__ int lin_right_locate_3d(double *data, int size, double to_locate, int num_adsorbates, int adsorbate, int dim);

__device__ __host__ int left_locate_3d(double *data, int lbound, int rbound, double to_locate, int num_adsorbates, int adsorbate, int dim, int depth);
__device__ __host__ int right_locate_3d(double *data, int lbound, int rbound, double to_locate, int num_adsorbates, int adsorbate, int dim);

#endif
