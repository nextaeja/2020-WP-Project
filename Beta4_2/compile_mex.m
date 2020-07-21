function  compile_mex()
%COMPILE_MEX Compile all the mex functions

mexcuda Setup/allocate_all_arrays.cu MEX_helpers/complex.cu
mexcuda MEX_helpers/print_CUDA_array.cu MEX_helpers/complex.cu
mexcuda MEX_helpers/print_complex_CUDA_array.cu MEX_helpers/complex.cu
mexcuda PropagationAlgorithms/compute_expk.cu MEX_helpers/complex.cu

%mexcuda Setup/print_complex_array.cu Setup/helper.cu MEX_helpers/complex.cu;
%mexcuda Setup/free_array.cu Setup/helper.cu;
%mexcuda PropagationAlgorithms/potential_propagation.cu Setup/helper.cu;
%mexcuda Setup/cuda_setup_dynamic_potential_2.cu Setup/helper.cu;
%mex Setup/setup_dynamic_potential.c;
%mexcuda Setup/print_potential.cu Setup/helper.cu;

end

