% Split Operator - O(dt^3) - V split
% TODO: Update for time varying V
function psiStepped = SplitOperatorStep_exp_3rdOrder_VSplit()
    global psi V mass hBar kSquared dt;
    
    expV = exp((-1i*(dt/2)/hBar)*V);
    expK = exp((-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass)));
    
    psiVStepHalf = expV.*psi;
    
    psiVStepHalfFT = fftn(psiVStepHalf);
    psiKStepFT = expK.*psiVStepHalfFT;
    psiKStep = ifftn(psiKStepFT);
    
    psiVStep = expV.*psiKStep;
    
    psiStepped = psiVStep;
end