#ifndef __COMPLEX__
#define __COMPLEX__

#include <math.h>

// Structure for a complex number
typedef struct complex {
	double real;
	double imag;

	__device__ __host__ complex(void): real(0), imag(0) {}
	__device__ __host__ complex(double x, double y): real(x), imag(y) {}

	// Multiplication
	__device__ __host__ complex operator* (const complex& z) {
		return complex(real*z.real - imag*z.imag, imag*z.real + real*z.imag);
	}
	__device__ __host__ complex operator* (const double& x) {
		return complex(real*x, imag*x);
	}

	// Addition
	__device__ __host__ complex operator+ (const complex& z) {
		return complex(real+z.real, imag+z.imag);
	}

	// Subtraction
	__device__ __host__ complex operator- (const complex& z) {
		return complex(real-z.real, imag-z.imag);
	}

	// Unary minus
	__device__ __host__ complex operator- () {
		return complex(-real, -imag);
	}

	// Compute exponential as exp(z) = exp(x+i*y) = exp(x)*(cos(y)+i*sin(y))
	__device__ __host__ complex _cexp(void) {
		return complex(cos(imag), sin(imag)) * exp(real);
	}
} complex;

const complex I = complex(0, 1);

__global__ void cexp(complex *zs, size_t size);

#endif
