#include "cuda_helper.h"

void _handle_cuda_error(cudaError_t code, const char *file, int line) {
	if (code != cudaSuccess) {
		char err_msg[500];
		sprintf(err_msg, "Error '%s' occurred in file '%s'@%d\n", cudaGetErrorString(code), file, line);
		mexErrMsgIdAndTxt("SplitOperator:CUDA:FFT", err_msg);
	}
}

void _handle_cudafft_error(cufftResult code, const char *file, int line) {
	if (code != CUFFT_SUCCESS) {
		char err_msg[500];
		sprintf(err_msg, "Cuda FFT error occurred in file '%s'@%d\n", file, line);
		mexErrMsgIdAndTxt("SplitOperator:CUDA:FFT", err_msg);
	}
}
