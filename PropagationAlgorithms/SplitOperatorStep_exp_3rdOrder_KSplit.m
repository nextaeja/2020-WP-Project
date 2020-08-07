% Split Operator - O(dt^3) - K split
function psiStepped = SplitOperatorStep_exp_3rdOrder_KSplit()
    global psi V mass hBar kSquared dt;
    
    expV = exp((-1i*dt/hBar)*V);
    expK = exp((-1i*(dt/2)/hBar)*(-hBar^2*-kSquared/(2*mass)));
    
    psiFT = fftn(psi);
    psiKStepHalfFT = expK.*psiFT;
    psiKStepHalf = ifftn(psiKStepHalfFT);
    
    psiVStep = expV.*psiKStepHalf;
    
    psiVStepFT = fftn(psiVStep);
    psiKStepFT = expK.*psiVStepFT;
    
    psiKStep = ifftn(psiKStepFT);
    
    psiStepped = psiKStep;
end