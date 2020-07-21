% Split Operator - O(dt^3) - V split
% TODO: Update for time varying V
% TODO: 

function psiStepped = SplitOperatorStep_exp_3rdOrder_VSplit_TimeDependent(potential,sp,wp,time)
   
    % Update potential V to V(t)
    UpdateBrownianMotionGaussians(sp, sp.decayType, sp.alpha, sp.xSigma, sp.ySigma, sp.gaussPeakVal, sp.wellDepth, time);
    
    expV = exp((-1i*(sp.dt/2)/SIUnits.hBar)*potential);
    expK = exp((-1i*sp.dt/SIUnits.hBar)*(-SIUnits.hBar^2*-sp.kSquared/(2*wp.mass)));
    
    psiVStepHalf = expV.*wp.psi;
    psiVStepHalfFT = fftn(psiVStepHalf);
    psiKStepFT = expK.*psiVStepHalfFT;
    psiKStep = ifftn(psiKStepFT);
    
    %%%%TODO:INTESTIGATE%SPLITTING%TIMES-V(dt)%vs%V(dt/2)%etc...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update potential V to V(t + dt) - see http://phys.au.dk/fileadmin/site_files/quscope/Haru_split_operator.pdf
    UpdateBrownianMotionGaussians(sp,sp.decayType, sp.alpha, sp.xSigma, sp.ySigma, sp.gaussPeakVal, sp.wellDepth, time + sp.dt);
    expV = exp((-1i*(sp.dt/2)/SIUnits.hBar)*potential);
    
    psiVStep = expV.*psiKStep;
    
    psiStepped = psiVStep;
end