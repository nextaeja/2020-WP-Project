/* Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
 * Copyright (c) 2018, Francis Haghighi-Daly 
 * All rights reserved.
 * This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.*/

#ifndef __COMPLEX__
#define __COMPLEX__

#include <math.h>

typedef double2 myComplex;

// Overloaded addition
__device__ myComplex element_complex_add(myComplex z1, myComplex z2);
__device__ myComplex element_complex_add(double x, myComplex z);
__device__ myComplex element_complex_add(myComplex z, double x);

// Overloaded multiplication
__device__ myComplex element_complex_mul(myComplex z1, myComplex z2);
__device__ myComplex element_complex_mul(double x, myComplex z);
__device__ myComplex element_complex_mul(myComplex z, double x);

__device__ myComplex element_complex_exp(myComplex z);

__global__ void complex_scale_shift(myComplex *zs, double scale, double shift, bool opposite, size_t size);
__global__ void complex_scale(myComplex *zs, double scale, size_t size);
__global__ void complex_exp(myComplex *zs, size_t size);
__global__ void complex_mul(myComplex *zs, myComplex *ws, size_t size);

#endif
