function UpdateBrownianMotionGaussians(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, t)
    global gaussianPositions it;
    
    

    x0=gaussianPositions(:,1,it+1);
    
    y0 = squeeze(gaussianPositions(:,2,it+1));
    
    SetupDynamicGaussianPotential_mex(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0)
end