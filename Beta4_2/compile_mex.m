function  compile_mex()
%COMPILE_MEX Compile all the mex functions

mexcuda Setup/allocate_all_arrays.cu MEX_helpers/complex.cu -lcuda
mexcuda Setup/cuda_setup_dynamic_potential_2.cu MEX_helpers/complex.cu -lcuda
mexcuda Setup/free_array.cu Setup/helper.cu -lcuda

mexcuda MEX_helpers/print_CUDA_array.cu MEX_helpers/complex.cu -lcuda
mexcuda MEX_helpers/print_complex_CUDA_array.cu MEX_helpers/complex.cu -lcuda
mexcuda MEX_helpers/print_complex_CUDA_array.cu MEX_helpers/complex.cu -lcuda
mexcuda MEX_helpers/copy_complex_array.cu MEX_helpers/complex.cu -lcuda
mexcuda MEX_helpers/cmp_complex_matlab_CUDA.cu MEX_helpers/complex.cu -lcuda

mexcuda PropagationAlgorithms/compute_expk.cu MEX_helpers/complex.cu -lcuda
mexcuda PropagationAlgorithms/compute_expv.cu MEX_helpers/complex.cu -lcuda
mexcuda PropagationAlgorithms/compute_fft_step.cu MEX_helpers/complex.cu -lcuda -lcufft
mexcuda PropagationAlgorithms/mex_split_operator_step_3rd_vsplit_time_dependent.cu MEX_helpers/complex.cu MEX_helpers/interpolation1d.cu Setup/cuda_setup_dynamic_potential.cu -lcuda -lcufft

end

