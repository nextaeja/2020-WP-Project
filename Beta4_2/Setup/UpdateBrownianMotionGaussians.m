function UpdateBrownianMotionGaussians(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, t)
    global gaussianPositions it;
    
    

    x0=gaussianPositions(:,1,it+1);
    
    y0 = squeeze(gaussianPositions(:,2,it+1));
    
    % Setup the potential in the "old" and "correct" way
    tic
    SetupDynamicGaussianPotential_mex(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0)
    standardTime = standardTime + toc;
    correctV = V;
    
    % Run the cuda function
    tic
    cuda_setup_dynamic_potential_2(V_ptr, z_offset_ptr, x0_ptr, y0_ptr, decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0, dx, dy, dz, A, nx, ny, nz);
    cudaTime = cudaTime + toc;
    
    nCalls = nCalls + 1;
end