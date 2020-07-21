function V = UpdateBrownianMotionGaussians(sp, decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, t)
    global ny gaussianPositions gaussianPositionsTimes numAdsorbates;
    
    tQuery = t;
    
    x0 = zeros(sp.numAdsorbates, 1);
    y0 = zeros(sp.numAdsorbates, 1);
    
    for adsorbateNum = 1:sp.numAdsorbates
        x0(adsorbateNum) = interp1(squeeze(gaussianPositionsTimes(adsorbateNum, 1,:)), squeeze(gaussianPositions(adsorbateNum, 1,:)), tQuery);	% x-centre of Gaussian

        if sp.ny ~= 1 % If 1D then interp1 will not interpolate because all values will be 0
            y0(adsorbateNum) = interp1(squeeze(gaussianPositionsTimes(adsorbateNum, 2,:)), squeeze(gaussianPositions(adsorbateNum, 2,:)), tQuery);	% y-centre of Gaussian
        else
            y0(adsorbateNum) = gaussianPositions(adsorbateNum, 2,1);
        end
    end
    
    V = SetupDynamicGaussianPotential(sp, decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0);
end