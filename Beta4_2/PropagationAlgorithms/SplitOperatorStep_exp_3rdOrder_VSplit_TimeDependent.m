% Split Operator - O(dt^3) - V split
% TODO: Update for time varying V
function psiStepped = SplitOperatorStep_exp_3rdOrder_VSplit_TimeDependent(t, V_ptr, z_offset_ptr, x0_ptr, y0_ptr, expV_ptr, expK_ptr, expK, psi_ptr)
    global psi V mass hBar kSquared dt decayType A eV ps nx ny nz;
   
    %decayType = 2; %%% 1 = exponential repulsive. 2 = Morse attractive. 3 = Morse-like (needs alpha parameter input too!)
    alpha = 2; % Only needed for Morse-like potential. alpha = 2 gives Morse potential. alpha = 0 gives exponential potential.
    xSigma = 3*(5.50/6)*A;       % x standard deviation
    ySigma = 3*(5.50/6)*A;       % y standard deviation
    gaussPeakVal = 3*1.61*A;   % peak value of Gaussian
    wellDepth = 10e-3*eV;
    
    % Update potential V to V(t)
    UpdateBrownianMotionGaussians(V_ptr, z_offset_ptr, x0_ptr, y0_ptr, decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth, t);
    
    expV = exp((-1i*(dt/2)/hBar)*V);
    compute_expv(V_ptr, hBar, dt, nx*ny*nz);
    cmp_complex_matlab_CUDA(gather(expV), V_ptr, 1e-4, nx*ny*nz);
    
    psiVStepHalf = expV.*psi;
    psiVStepHalfFT = fftn(psiVStepHalf);
    psiKStepFT = expK.*psiVStepHalfFT;
    psiKStep = ifftn(psiKStepFT);
    
    % Reverse the x-y
    bangalla = permute(psiVStepHalf, [3 1 2]);
    yadda = fftn(bangalla);
    
    compute_fft_step(V_ptr, expK_ptr, psi_ptr, nx, ny, nz);
    print_complex_CUDA_array(psi_ptr, nx, ny, nz);
    cmp_complex_matlab_CUDA(gather(yadda), psi_ptr, 1e-5, nx*ny*nz, 1);
    error('test fft');
    
    %%%%TODO:INTESTIGATE%SPLITTING%TIMES-V(dt)%vs%V(dt/2)%etc...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update potential V to V(t + dt) - see http://phys.au.dk/fileadmin/site_files/quscope/Haru_split_operator.pdf
    UpdateBrownianMotionGaussians(V_ptr, z_offset_ptr, x0_ptr, y0_ptr, decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth, t + dt);
    expV = exp((-1i*(dt/2)/hBar)*V);
    
    psiVStep = expV.*psiKStep;
    
    psiStepped = psiVStep;
    error('stop');
end