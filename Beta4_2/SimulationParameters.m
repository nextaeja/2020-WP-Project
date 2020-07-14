function sp = SimulationParameters()
        sp.decayType = 3;
        sp.alpha = 2; % Only needed for Morse-like potential. alpha = 2 gives Morse potential. alpha = 0 gives exponential potential.
        sp.propagationMethod = 5;
        sp.numAdsorbates = 30;

        si = SIUnits;
        sp.lx = 90*si.A;
        sp.ly = 90*si.A;
        sp.lz = 90*si.A;
        % Acceptable error in wavepacket norm
        sp.eps = 1e-6;
        
        % Time properties
        sp.dt0 = 0.01*si.ps;      % Initial timestepsp. Units = s
        sp.tStart = 0*si.ps;      % Units = s
        sp.tFinish = 12*si.ps+sp.dt0;    % Units = s

        sp.savingSimulationRunning = false;
        sp.savingSimulationEnd = false;
        sp.realTimePlotting = true;
        sp.displayAdsorbateAnimation = false;

        sp.numPsiToSave = 5;
        sp.numGfxToSave = 20;

        sp.nx = 64;
        sp.ny = 64;
        sp.nz = 256;
        
        % Brownian motion
        sp.xSigma = 3*(5.50/6)*si.A;       % x standard deviation
        sp.ySigma = 3*(5.50/6)*si.A;       % y standard deviation
        sp.gaussPeakVal = 3*1.61*si.A;   % peak value of Gaussian
        sp.wellDepth = 10e-3*si.eV;

end