% Ocean Haghighi-Daly - Part III project
% Working attempt at Split-Operator method of wavepacket propagation
% Done in 3D
%
% Some code taken from Andrew Alderwick

% Beta 4_x - see 'details' file for version details

% Units: SI Units:
% hBar = 1.054571800e-34 Js. c = 299792458 m/s. eV = 1.6021766208e-19 J.
% 1amu =  1.660539040e-27 kg. 1Å = 1e-10m. 1ps = 1e-12s.
% 3He mass = 3.0160293amu
% 1Å ~ Atomic diameter
                                                                                                                                                                                                                                                                                                                                                                                                                                           
function WavepacketPropagation_beta4_2
    tbegin = tic;  
%===START=TIMING========================================================================================%
    tic
    
%===ADD=PATHS=TO=OTHER=FUNCTIONS========================================================================%
    % genpath gets paths to all folders in this folder
    addpath(genpath(pwd));
    
%===SET=UP=VARIABLES====================================================================================%
    global A ps eV hBar kSquared mass nx ny nz dx dy dz lx ly lz eps tStart tFinish notifySteps gfxSteps psi psi0 dt0 dt savingSimulationRunning savingDirectory propagationMethod numAdsorbates decayType custpot zOffset pathfile Browniefile savingBrownianPaths it numIterations gaussianPositions
    
    gpuDevice(1);
    SetupSIUnits();
    
    % Setup lengths. Units = m
    lx = 90*A;
    ly = 90*A;
    lz = 90*A;
    
    % Setup grid - use powers of 2 for quickest FFT
    nx = 128;
    ny = 128;
    nz = 512;
    
    % Acceptable error in wavepacket norm
    eps = 1e-4;
    
    % Time properties
    dt0 = 0.01*ps;      % Initial timestep. Units = s
    tStart = 0*ps;      % Units = s
    tFinish = 12*ps;    % Units = s
        
    realTimePlotting = true; %also determines if graphics are saved for SavingSimulationRunning
    savingSimulationRunning = false;
    savingSimulationEnd = false;
    displayAdsorbateAnimation = true;
    savingBrownianPaths=false;
    Browniefile="PotentialFiles/brownianpaths.txt";
    
    saveLocation =""; %if "" it will save in beta 4_2
    
    numPsiToSave = 1;%not used for prop method 6
    numGfxToSave = 20; %for prop 6, psi saved at the same time, psi can be saved while disabling graphics by disabling realTimePlotting 
    
    % Propagation method: 1 = RK4Step. 2 = Split Operator O(dt^2). 3 = Split Operator O(dt^3), K split. 4 = Sp. Op. O(dt^3), V split. 5 = Sp.Op. O(dt^3), V
    % split, time dependent. 6 = MEXCUDA while loop with Sp.Op. O(dt^3), V split, time dependent.
    % 7 = matlab and MEX Sp.Op. O(dt^3), V split, time dependent. 8 = MEX Sp. Op. O(dt^3) V split,
    % time dependent. 
    propagationMethod = 6;
    numAdsorbates = 30;
    
    custompaths = true;
    pathfile = "PotentialFiles/brownianpaths.txt";
        
    if (custompaths == true) %autoredetermines what num adsorbates is for a custom potential as different numAdsorbates between customPotential and program result in strange errors (perhaps disable this section for 2D scattering as it has not yet been tested)
        fid = fopen(pathfile,'rt');
        tLines = fgets(fid);
        numCols = numel(strfind(tLines,' '))/2;
        fclose(fid);
        numAdsorbates =  numCols;
    end
    
    
    numSteps = round(tFinish/dt0);
    
    
    decayType = 3; % 1 = exponential repulsive. 2 = Morse attractive. 3 = Morse-like (needs alpha parameter input too!). 4=custom
    
    potfile = "PotentialFiles/potential.txt"; %for 4, text file containing floats for real and imaginary part of potential, seperated by lz. potential should be high to prevent tunelling over cyclic boundary  
    
    zOffset = -5*A; % Shift entire V away from boundary to stop Q.Tunneling through V %%% or to make custom potential go all the way to the surface
    
    alpha = 2; % Only needed for Morse-like potential. alpha = 2 gives Morse potential. alpha = 0 gives exponential potential.
    xSigma = 3*(5.50/6)*A;       % x standard deviation
    ySigma = 3*(5.50/6)*A;       % y standard deviation
    gaussPeakVal = 3*1.61*A;   % peak value of Gaussian
    wellDepth = 10e-3*eV;
    
    
     %End of inputs==========================================================%
     
    
    notifySteps = floor(numSteps/numGfxToSave);   % TODO: Change to notifytime. # steps after which to notify user of progress
    psiSaveSteps = floor(numSteps/numPsiToSave);
    
    if realTimePlotting&&(numGfxToSave ~= 0)
        gfxSteps = floor(numSteps/numGfxToSave);      % TODO: Change to gfxtime # steps after which, update graphics
    else
        if(propagationMethod == 6)
            if(numGfxToSave ~= 0)
                gfxSteps = floor(numSteps/numGfxToSave);
            else
                gfxSteps = numSteps;
            end
        else
            gfxSteps = 0;
        end
    end
    
    
    % Use integers for loop equality test. CARE: round will give you closest # to tFinish/dt and might be floor or ceiling value
    numIterations = round(tFinish/dt0);
    
    % Allocate CUDA arrays
    CUDA_pointers = SetupVariables();
    z_offset_ptr = CUDA_pointers(1);
    x0_ptr = CUDA_pointers(2);
    y0_ptr = CUDA_pointers(3);
    k_squared_ptr = CUDA_pointers(4);
    exp_v_ptr = CUDA_pointers(5);
    exp_k_ptr = CUDA_pointers(6);
    psi_ptr = CUDA_pointers(7);
    gauss_position_ptr = CUDA_pointers(8);
    
    if(~custompaths)
        SetupBrownianMotionGaussians(displayAdsorbateAnimation, realTimePlotting);%%%NaN bug caused by something in here %%% 
    else
       SetupSpecifiedPaths(displayAdsorbateAnimation, realTimePlotting);
    end
   
    %SetupZeroPotential();
    %SetupWedgePotential();
    
    if decayType==4
       f=fopen(potfile,'r');
       custpot = fscanf(f,'%f'); 
       fclose(f);
    end
    
    
    % it = iteration. Starts at 0, represents the iteration CURRENTLY being carried out.
    it = 0;
    
    UpdateBrownianMotionGaussians(decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth, tStart);
    %SetupStaticGaussianPotential(decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth);

    % Initialises psi0
    SetupInitialWavefunction();
    
    % Copy the initialized function into the allocated space
    % gather is necessary to pass a gpuArray into RAM
    copy_CUDA_complex_array(psi_ptr, gather(psi0), nx*ny*nz);
    copy_CUDA_array(gauss_position_ptr, gaussianPositions, numAdsorbates*2*numIterations);
    
    SetupGraphicsVariables();

%===SET=UP=COMPLETE=====================================================================================%
    fprintf('Variable setup complete. Time: %.6fs.\n', toc);
    
%===SAVING=RESULTS=====================================================================================%
    % Variable dictating whether simulation saved in 'Simulations' folder
    % or not 
    if(saveLocation=="")
        savingDirectory = strcat(pwd,'\SavedSimulation');
    else
        savingDirectory=saveLocation;
    end
    
    % Create unique simulation folder to store results
    if(savingSimulationRunning || savingSimulationEnd)
        CreateNewSimulationFolder();
    end
    
    % Save initialisation data
    if(savingSimulationRunning || savingSimulationEnd)
        SaveInitialisationData(alpha, xSigma, ySigma, gaussPeakVal, wellDepth, savingDirectory);
    end
    
    
%===RUN=SIMULATION======================================================================================%
    % Initial conditions
    psi = psi0;
    dt = dt0;
    t = tStart;
        
    % Compute the value of expK
    % This is constant unless the value of dt is changed
    expK = exp((-1i*dt/hBar)*(-hBar^2*-kSquared/(2*mass)));
    compute_expk(exp_k_ptr, k_squared_ptr, hBar, dt, mass, kSquared, nx*ny*nz);
    cmp_complex_matlab_CUDA(expK, exp_k_ptr, eps, nx, ny, nz);
    
    % Loop iteratively until tFinish reached
    standardTime = 0.0;
    cudaTime = 0.0;
    nCalls = 0;
    t = tStart;
    
    if(propagationMethod==6)
        for i=1:numGfxToSave
            if(realTimePlotting)
                UpdateGraphics(t, it)  
                if savingSimulationRunning
                    SaveSimulationRunning(t);
                end
            end
            mexcudawhile(exp_v_ptr, z_offset_ptr, gauss_position_ptr, x0_ptr, y0_ptr, exp_k_ptr, psi_ptr, nx, ny, nz, decayType, A, eV, hBar, dt, dx, dy, dz, gfxSteps,t,alpha,it,numAdsorbates);         
            it = (i)*(gfxSteps);
            t=it*dt+tStart;
            
            
            psi = copy_from_CUDA_complex(psi_ptr, nx, ny, nz);
            if savingSimulationRunning
                saveName = strcat('psi_t', num2str(t/ps), '.mat');
                save(saveName, 'psi');
            end
            
            % Should move this line to the cuda while ?=================%
            totProb = sum(sum(sum(psi.*conj(psi))));
        
            % Check unitarity (Note that the mex while loop only checks for unitarity when it makes a graphic, and not every iteration like the others)
            if abs(totProb - 1) > eps
                 fprintf(1, 'Step %d incomplete - unitarity error caused by previous step: %d. Time (%.3f ps, %.3f s): (unitarity %.7f)\n', it, it - 1, (it - 1)*dt/ps, toc, totProb);
                 error("unitary error") %now terminates on this rather than restarting
            end
        end
        %it = it + 1; %to retain consistency for final graphic
    else

    while(it < numIterations)  
        % Total probability
        totProb = sum(sum(sum(psi.*conj(psi))));
        
        % Check unitarity
        if abs(totProb - 1) > eps
            fprintf(1, 'Step %d incomplete - unitarity error caused by previous step: %d. Time (%.3f ps, %.3f s): (unitarity %.7f)\n', it, it - 1, (it - 1)*dt/ps, toc, totProb);
            error("unitary error") %now terminates on this rather than restarting
        else
            % Notify user if necessary. it - 1 as step is NOT complete yet. it - 1 is complete.
            if notifySteps > 0 && mod(it, notifySteps)== 0
                fprintf(1, 'Step %d complete (%.3f ps, %.3f s): propagate wpkt (unitarity %.7f)\n', it, t/ps, toc, totProb);
                fprintf("MATLAB time %.3f, CUDA time %.3f, speedup x%.3f\n", standardTime, cudaTime, standardTime / cudaTime);
            end
            
            % Produce graphics if asked and if correct # of steps has passed
            if gfxSteps > 0 && mod(it, gfxSteps) == 0p
                UpdateGraphics(t, it)
                
                if savingSimulationRunning
                    SaveSimulationRunning(t);
                end
            end
            % Save psi if necessary
            if savingSimulationRunning && psiSaveSteps > 0 && mod(it - 1, psiSaveSteps) == 0
                saveName = strcat('psi_t', num2str(t/ps), '.mat');
                save(saveName, 'psi');
            end
            
%===========STEP=FORWARD================================================================================%
            
            switch propagationMethod
                case 1
                    psi = RK4Step();
                case 2
                    psi = SplitOperatorStep_exp();
                case 3
                    psi = SplitOperatorStep_exp_3rdOrder_KSplit();
                case 4
                    psi = SplitOperatorStep_exp_3rdOrder_VSplit();
                case 5 %noncuda original iteration method
                    psi = SplitOperatorStep_exp_3rdOrder_VSplit_TimeDependent(t, expK);
                case 7
                    tic;
                    psi = SplitOperatorStep_exp_3rdOrder_VSplit_TimeDependent(t, expK);
                    standardTime = standardTime + toc;

                    tic;
                    mex_split_operator_step_3rd_vsplit_time_dependent(t, exp_v_ptr, z_offset_ptr, gauss_position_ptr, x0_ptr, y0_ptr, exp_k_ptr, psi_ptr, nx, ny, nz, decayType, A, eV, hBar, dt, dx, dy, dz, it, numAdsorbates); %currently always has alpha=2
                    mex_psi = copy_from_CUDA_complex(psi_ptr, nx, ny, nz);
                    cudaTime = cudaTime + toc;
                    nCalls = nCalls + 1;
                case 8 %Cuda only (which doesn't seem to exist..?)
                    tic;
                    mex_split_operator_step_3rd_vsplit_time_dependent(t, exp_v_ptr, z_offset_ptr, gauss_position_ptr, x0_ptr, y0_ptr, exp_k_ptr, psi_ptr, nx, ny, nz, decayType, A, eV, hBar, dt, dx, dy, dz, it, numAdsorbates); %currently always has alpha=2
                    psi = copy_from_CUDA_complex(psi_ptr, nx, ny, nz);
                    cudaTime = cudaTime + toc;
                    nCalls = nCalls + 1;
            end
            % Iteration it complete. t is now t + dt
            t = t + dt;
            
            % Iterate counters at end. it complete, next step is it + 1
            it = it + 1;
        end
    end %While
    end
    
    if(propagationMethod == 7)
        psiCudaErrors = (psi - mex_psi);
        psiPlaceholder = psiCudaErrors./mex_psi;
        psiErrorSize = size(psiPlaceholder);
        
        psiErrorTotal = sum(psiPlaceholder, 'all');
        psiErrorTotal = abs(psiErrorTotal);
        
        psiErrorSize = prod(psiErrorSize);
        
        avgPsiError = psiErrorTotal/psiErrorSize;
        avgPsiError = 100 * avgPsiError;
    end
    
    tEnd = toc(tbegin);  
    %}
    % Tell user run is complete
    % Note, finalIteration = it - 1 as it starts counting at 1. t starts at 0 though, and t represents the time just after the last iteration, so tFinal = t
    fprintf('Run Complete.\nNumber of iterations = %d\nFinal simulation time = %.16e\n', it , t);
    if propagationMethod == 7	
        fprintf("MATLAB time %.3f, CUDA time %.3f, speedup x%.3f, %% error in CUDA vs Matlab %.5f\n", standardTime, cudaTime, standardTime / cudaTime, avgPsiError);
    else
        fprintf('Total Run Time %.3f\n', tEnd);
    end
    % Force graphics update so psi_final is displayed
    if (gfxSteps > 0 &&realTimePlotting)
        UpdateGraphics(t, it)
    end
    
    % If saving simulation, print end data to file
    if(savingSimulationRunning || savingSimulationEnd)
       SaveSimulationEndData(t, it);
    end
    
    %{
    exp_v_ptr
    z_offset_ptr
    gauss_position_ptr
    x0_ptr
    y0_ptr
    exp_k_ptr
    psi_ptr
    %}
    
    % Free previously allocated memory in MEX files
    free_array(0, size(CUDA_pointers, 2), [], CUDA_pointers);
end
