function UpdateBrownianMotionGaussians(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, t)
    global gaussianPositions it;
    
    

    x0=gaussianPositions(:,1,it+1);
    
    y0 = squeeze(gaussianPositions(:,2,it+1));
    
    % Setup the potential in the "old" and "correct" way
    tic
    SetupDynamicGaussianPotential_mex(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0)
    standardTime = standardTime + toc;
    correctV = V;
    
    % Run the rewritten function
    tic
    c_potential = setup_dynamic_potential(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0, dx, dy, dz, A, nx, ny, nz);
    cTime = cTime + toc;
    
    % Run the cuda function
    tic
    cuda_potential = cuda_setup_dynamic_potential(decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0, dx, dy, dz, A, nx, ny, nz);
    cudaTime = cudaTime + toc;
    
    % Check for correctness
    assert(isequal(size(correctV), size(c_potential)));
    assert(isequal(correctV, c_potential));
    assert(isequal(size(correctV), size(cuda_potential)));
    assert(isequal(correctV, cuda_potential));
    
    nCalls = nCalls + 1;
end