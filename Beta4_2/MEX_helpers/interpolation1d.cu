#include "interpolation1d.h"

// Save in dev_x0 and dev_y0 the position value at time t_query with linear interpolation
__global__ void interpolate1d_adsorbate_positions(double *gaussian_times, double *gaussian_positions, int num_adsorbates, int num_gaussian_positions, double *dev_x0, double *dev_y0, double t_query, int ny) {
	int adsorbate = blockIdx.x * blockDim.x + threadIdx.x;

	while (adsorbate < num_adsorbates) {
		// Get bounds for time interval around t_query
		int min_time_idx = 0;//left_locate_3d(gaussian_times, 0, num_gaussian_positions-1, t_query, num_adsorbates, adsorbate, 0);
		int max_time_idx = 0;//right_locate_3d(gaussian_times, 0, num_gaussian_positions-1, t_query, num_adsorbates, adsorbate, 0);

		// Handle case of a out of range value
		if (min_time_idx < 0 || max_time_idx < 0) {
			dev_x0[adsorbate] = 0.0;
			dev_y0[adsorbate] = 0.0;
			break;
		}

		// Extract the values for the given adsorbate
		double min_time = get_gaussian_adsorbate(gaussian_times, min_time_idx, 0, adsorbate, num_adsorbates);
		double max_time = get_gaussian_adsorbate(gaussian_times, max_time_idx, 0, adsorbate, num_adsorbates);
		double min_xpos = get_gaussian_adsorbate(gaussian_positions, min_time_idx, 0, adsorbate, num_adsorbates);
		double max_xpos = get_gaussian_adsorbate(gaussian_positions, max_time_idx, 0, adsorbate, num_adsorbates);
		double min_ypos = get_gaussian_adsorbate(gaussian_positions, min_time_idx, 1, adsorbate, num_adsorbates);
		double max_ypos = get_gaussian_adsorbate(gaussian_positions, max_time_idx, 1, adsorbate, num_adsorbates);

		// If the time was tabulated return the result
		if (min_time_idx == max_time_idx) {
			dev_x0[adsorbate] = min_xpos;
			dev_y0[adsorbate] = min_ypos;
		} else {
			// Interpolate x position
			double gradient = (max_xpos - min_xpos) / (max_time - min_time);
			dev_x0[adsorbate] = gradient * (t_query - min_time) + min_xpos;

			// Interpolate y position
			if (ny == 1) { // Handle 1D case
				dev_y0[adsorbate] = get_gaussian_adsorbate(gaussian_positions, 0, 1, adsorbate, num_adsorbates);
			} else {
				double gradient = (max_ypos - min_ypos) / (max_time - min_time);
				dev_y0[adsorbate] = gradient * (t_query - min_time) + min_ypos;
			}
		}

		adsorbate += blockDim.x * gridDim.x;
	}
}

// The gaussian position array is a 3D one, return the correct value for a given
//	adsorbate number and dimension (x or y)
__device__ __host__ double get_gaussian_adsorbate(double *data, int idx, int dim, int adsorbate, int num_adsorbates) {
	if (idx < 0) {
		return -1.0;
	}

	int tot_idx = adsorbate + dim*num_adsorbates + idx*num_adsorbates*NUM_GAUSSIAN_ADSORBATE_DIMENSIONS;

	return data[tot_idx];
}

// With data being a sorted array, return the largest element smaller than to_locate
__device__ __host__ int left_locate_3d(double *data, int lbound, int rbound, double to_locate, int num_adsorbates, int adsorbate, int dim, int depth) {
	if (depth > 100) {
		mexErrMsgIdAndTxt("SplitOperator:Interpolate:LeftLocate", "Hit max recursion depth in left locate\n");
		mexPrintf("Bangalla");
		return -1;
	} else {
		depth++;
	}

	if (rbound - lbound == 1) return lbound;

	// Find mid point
	int mid_point = lbound + (rbound - lbound) / 2;

	// Extract the values
	double lvalue = get_gaussian_adsorbate(data, lbound, dim, adsorbate, num_adsorbates);
	double rvalue = get_gaussian_adsorbate(data, rbound, dim, adsorbate, num_adsorbates);
	double mid_value = get_gaussian_adsorbate(data, mid_point, dim, adsorbate, num_adsorbates);
	//mexPrintf("%e) %e %e %e\n", to_locate, lvalue, mid_value, rvalue);


	// Handle the case of out of bound value
	if (lvalue > to_locate || rvalue < to_locate) return -1;

	// Handle value found at extremities
	if (lvalue == to_locate) return lbound;
	if (rvalue == to_locate) return rbound;

	if (mid_value == to_locate) {	// Value found in the middle
		mexPrintf("Tho bekkato\n");
		return mid_point;
	} else if (to_locate < mid_value) {
		//mexPrintf("Going left\n");
		return left_locate_3d(data, lbound, mid_point, to_locate, num_adsorbates, adsorbate, dim, depth);
	} else {
		//mexPrintf("Going right\n");
		return left_locate_3d(data, mid_point+1, rbound, to_locate, num_adsorbates, adsorbate, dim, depth);
	}
}

// With data being a sorted array, return the smallest element larger than to_locate
__device__ __host__ int right_locate_3d(double *data, int lbound, int rbound, double to_locate, int num_adsorbates, int adsorbate, int dim) {
	if (rbound - lbound == 1) return rbound;
	static int depth = 0;
	if (depth > 10000) {
		mexPrintf("Hit max recursion depth in right locate\n");
		return -1;
	} else {
		depth++;
	}

	// Find mid point
	int mid_point = lbound + (rbound - lbound) / 2;

	// Extract the values
	double lvalue = get_gaussian_adsorbate(data, lbound, dim, adsorbate, num_adsorbates);
	double rvalue = get_gaussian_adsorbate(data, rbound, dim, adsorbate, num_adsorbates);
	double mid_value = get_gaussian_adsorbate(data, mid_point, dim, adsorbate, num_adsorbates);

	// Handle the case of out of bound value
	if (lvalue > to_locate || rvalue < to_locate) return -1;

	// Handle value found at extremities
	if (lvalue == to_locate) return lbound;
	if (rvalue == to_locate) return rbound;

	if (mid_value == to_locate) {	// Value found in the middle
		return mid_point;
	} else if (to_locate < mid_value) {
		return right_locate_3d(data, lbound, mid_point, to_locate, num_adsorbates, adsorbate, dim);
	} else {
		return right_locate_3d(data, mid_point+1, rbound, to_locate, num_adsorbates, adsorbate, dim);
	}
}
