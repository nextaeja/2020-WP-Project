#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"
#include "../Setup/cuda_setup_dynamic_potential.h"


const int NUM_GAUSSIAN_ADSORBATE_DIMENSIONS = 2;
