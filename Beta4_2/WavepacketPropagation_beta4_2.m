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

%===START=TIMING========================================================================================%
    tic
    
%===ADD=PATHS=TO=OTHER=FUNCTIONS========================================================================%
    % genpath gets paths to all folders in this folder
    addpath(genpath(pwd));
    
%===SET=UP=VARIABLES====================================================================================%
    global A ps eV nx ny nz lx ly lz eps tStart tFinish notifySteps gfxSteps psi psi0 dt0 dt savingSimulationRunning savingDirectory propagationMethod numAdsorbates decayType
    
    SetupSIUnits();
    
    % Setup lengths. Units = m
    lx = 90*A;
    ly = 90*A;
    lz = 90*A;
    
    % Setup grid - use powers of 2 for quickest FFT
    nx = 64;
    ny = 64;
    nz = 256;
    
    % Acceptable error in wavepacket norm
    eps = 1e-6;
    
    % Time properties
    dt0 = 0.01*ps;      % Initial timestep. Units = s
    tStart = 0*ps;      % Units = s
    tFinish = 12*ps;    % Units = s
        
    savingSimulationRunning = false;
    savingSimulationEnd = true;
    realTimePlotting = true;
    displayAdsorbateAnimation = false;

    numPsiToSave = 5;
    numGfxToSave = 20;
    numSteps = round(tFinish/dt0);
    
    notifySteps = floor(numSteps/numGfxToSave);   % TODO: Change to notifytime. # steps after which to notify user of progress
    psiSaveSteps = floor(numSteps/numPsiToSave);
    
    if realTimePlotting
        gfxSteps = floor(numSteps/numGfxToSave);      % TODO: Change to gfxtime # steps after which, update graphics
    else
        gfxSteps = 0;
    end
    
    % Propagation method: 1 = RK4Step. 2 = Split Operator O(dt^2). 3 = Split Operator O(dt^3), K split. 4 = Sp. Op. O(dt^3), V split. 5 = Sp.Op. O(dt^3), V
    % split, time dependent.
    propagationMethod = 5;
    
    SetupVariables();
    

    numAdsorbates = 30;
    SetupBrownianMotionGaussians(displayAdsorbateAnimation, realTimePlotting);%%%NaN bug caused by something in here
    
    tFinish=tFinish+dt; %%%NaN bug
   
    %SetupZeroPotential();
    %SetupWedgePotential();
    decayType = 2; % 1 = exponential repulsive. 2 = Morse attractive. 3 = Morse-like (needs alpha parameter input too!)
    alpha = 2; % Only needed for Morse-like potential. alpha = 2 gives Morse potential. alpha = 0 gives exponential potential.
    xSigma = 3*(5.50/6)*A;       % x standard deviation
    ySigma = 3*(5.50/6)*A;       % y standard deviation
    gaussPeakVal = 3*1.61*A;   % peak value of Gaussian
    wellDepth = 10e-3*eV;
    UpdateBrownianMotionGaussians(decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth, tStart);
    %SetupStaticGaussianPotential(decayType, alpha, xSigma, ySigma, gaussPeakVal, wellDepth);

    % Initialises psi0
    SetupInitialWavefunction();
    
    SetupGraphicsVariables();

%===SET=UP=COMPLETE=====================================================================================%
    fprintf('Variable setup complete. Time: %.6fs.\n', toc);
    
%===SAVING=RESULTS=====================================================================================%
    % Variable dictating whether simulation saved in 'Simulations' folder
    % or not 
    savingDirectory = strcat(pwd,'\SavedSimulation');
    
    % Create unique simulation folder to store results
    if(savingSimulationRunning || savingSimulationEnd)
        CreateNewSimulationFolder();
    end
    
    % Save initialisation data
    if(savingSimulationRunning || savingSimulationEnd)
        SaveInitialisationData(alpha, xSigma, ySigma, gaussPeakVal, wellDepth, savingDirectory);
    end
    
    % Make display full screen if saving simulation and gfxSteps > 0
    if(savingSimulationRunning && gfxSteps > 0)
        figure('units','normalized','outerposition',[0 0 1 1])
    end
    
%===RUN=SIMULATION======================================================================================%
    % Initial conditions
    psi = psi0;
    dt = dt0;
    t = tStart;
   
    
    % it = iteration. Starts at 1, represents the iteration CURRENTLY being carried out.
    it = 1;
    
    % Use integers for loop equality test. CARE: round will give you closest # to tFinish/dt and might be floor or ceiling value
    numIterations = round(tFinish/dt)-1;
    
    
    % Loop iteratively until tFinish reached
    while(it <= numIterations)  
        % Total probability
        totProb = sum(sum(sum(psi.*conj(psi))));
        
        % Check unitarity
        if abs(totProb - 1) > eps
            fprintf(1, 'Step %d incomplete - unitarity error caused by previous step: %d. Time (%.3f ps, %.3f s): (unitarity %.7f)\n', it, it - 1, (it - 1)*dt/ps, toc, totProb);
            disp(['Halving timestep to dt = ' num2str(dt/2) ' UNITS and restarting.'])
            errorMessage = strcat('Unitarity error after completing ', num2str(it - 1), 'step(s) - halving timestep. New dt = ', num2str(dt/2, '%.16e'));
            PrintErrorMessage(errorMessage, false); % Don't print in terminal - already done above.
            dt = dt/2;
            psi = psi0;
            it = 1;
            t = tStart;
            numIterations = round(tFinish/dt);
            notifySteps = floor(numIterations/numGfxToSave);
            gfxSteps = floor(numIterations/numGfxToSave);
            psiSaveSteps = floor(numIterations/numPsiToSave);
            
            % All items reset. Restarting simulation.
        else
            % Notify user if necessary. it - 1 as step is NOT complete yet. it - 1 is complete.
            if notifySteps > 0 && mod(it - 1, notifySteps)== 0
                fprintf(1, 'Step %d complete (%.3f ps, %.3f s): propagate wpkt (unitarity %.7f)\n', it - 1, t/ps, toc, totProb);
            end
            
            % Produce graphics if asked and if correct # of steps has passed
            if gfxSteps > 0 && mod(it - 1, gfxSteps) == 0
                UpdateGraphics(t, it - 1)
                
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
                case 5
                    psi = SplitOperatorStep_exp_3rdOrder_VSplit_TimeDependent(t);
            end
            % Iteration it complete. t is now t + dt
            t = t + dt;
            
            % Iterate counters at end. it complete, next step is it + 1
            it = it + 1;
        end
        
    end %While
    %numIterations=numIterations-1;
    % Tell user run is complete
    % Note, finalIteration = it - 1 as it starts counting at 1. t starts at 0 though, and t represents the time just after the last iteration, so tFinal = t
    fprintf('Run Complete.\nNumber of iterations = %d\nFinal simulation time = %.16e\n', it - 1, t);
    
    % Force graphics update so psi_final is displayed
    if gfxSteps > 0
        UpdateGraphics(t, it - 1)
    end
    
    % If saving simulation, print end data to file
    if(savingSimulationRunning || savingSimulationEnd)
       SaveSimulationEndData(t, it - 1);
    end
end
