% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function UpdateBrownianMotionGaussians(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, t)
    global gaussianPositions it;
    
    x0=gaussianPositions(:,1,it+1);
    
    y0 = squeeze(gaussianPositions(:,2,it+1));
    
    % Setup the potential in the "old" and "correct" way
    SetupDynamicGaussianPotential(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0)
end