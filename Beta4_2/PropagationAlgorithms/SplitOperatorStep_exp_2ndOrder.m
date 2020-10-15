% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

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