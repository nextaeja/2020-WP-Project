/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#ifndef __SETUP_DYNAMIC_POTENTIAL__
#define __SETUP_DYNAMIC_POTENTIAL__

#include "../MEX_helpers/complex.h"

void setup_dynamic_gaussian_potential(myComplex *dev_potential, double *dev_z_offset, double *dev_x0, double *dev_y0,
        int adsorbate_num, int nx, int ny, int nz, int decayType, int inParameterIfNeeded, double eV, double A,
        double dx, double dy, double dz);

#endif
