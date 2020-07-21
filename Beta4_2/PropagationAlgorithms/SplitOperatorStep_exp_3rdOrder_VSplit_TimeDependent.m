% Split Operator - O(dt^3) - V split
% TODO: Update for time varying V
function psiStepped = SplitOperatorStep_exp_3rdOrder_VSplit_TimeDependent(t, V_ptr, z_offset_ptr, x0_ptr, y0_ptr, expK)
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
    
    psiVStepHalf = expV.*psi;
    psiVStepHalfFT = fftn(psiVStepHalf);
    psiKStepFT = expK.*psiVStepHalfFT;
    psiKStep = ifftn(psiKStepFT);
    
    %%%%TODO:INTESTIGATE%SPLITTING%TIMES-V(dt)%vs%V(dt/2)%etc...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update potential V to V(t + dt) - see http://phys.au.dk/fileadmin/site_files/quscope/Haru_split_operator.pdf
    UpdateBrownianMotionGaussians(V_ptr, z_offset_ptr, x0_ptr, y0_ptr, decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth, t + dt);
    expV = exp((-1i*(dt/2)/hBar)*V);
    
    psiVStep = expV.*psiKStep;
    
    psiStepped = psiVStep;
    error('stop');
end