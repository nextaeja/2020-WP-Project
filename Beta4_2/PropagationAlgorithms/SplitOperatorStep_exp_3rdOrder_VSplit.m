% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

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