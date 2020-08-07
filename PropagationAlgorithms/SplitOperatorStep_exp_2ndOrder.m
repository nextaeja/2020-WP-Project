% Split Operator - O(dt^2) - equal K and V steps
function psiStepped = SplitOperatorStep_exp_2ndOrder()
    global psi V mass hBar kSquared dt;
    
    % -hBar^2/2m etc. MUST be INSIDE exponentials
    expV = exp((-1i*dt/hBar)*V);
    expK = exp((-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass)));
    
    psiVStep = expV.*psi;
    
    psiVStepFT = fftn(psiVStep);
    psiKStepFT = expK.*psiVStepFT;
    
    psiKStep = ifftn(psiKStepFT);
    
    psiStepped = psiKStep;
end