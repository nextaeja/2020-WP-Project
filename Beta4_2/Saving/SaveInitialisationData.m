% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function SaveInitialisationData(alphaIn, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, saveLocation)
    global startingDirectory serialNumber V; % Needed in function
    global nx ny nz lx ly lz dx dy dz tStart tFinish dt0 hBar c eV amu A ps propagationMethod;

    % Create file to save initialisation data in. 'w+' argument says open or create, and wipe any existing data...
    fileID = fopen('InitilisationData.txt', 'w+');

    fprintf(fileID, 'Simulation Serial Number: %05d\r\n\r\n', serialNumber); % need \r\n for new line...
    fprintf(fileID, 'Using program: %s\r\n', mfilename()); % mfilename gets name of currently running file
    fprintf(fileID, 'in folder: %s\r\n\r\n', startingDirectory);

    fprintf(fileID, 'Simulation Data:\r\n');
    fprintf(fileID, '----------------\r\n');

    fprintf(fileID, '(nx, ny, nz) = (%d, %d, %d)\r\n', nx, ny, nz);
    fprintf(fileID, '(lx, ly, lz) = (%g*A, %g*A, %g*A)\r\n', lx/A, ly/A, lz/A);
    fprintf(fileID, '(dx, dy, dz) = (%g*A, %g*A, %g*A)\r\n', dx/A, dy/A, dz/A);
    fprintf(fileID, 'Starting dt = dt0 = %gps\r\n', dt0/ps);
    fprintf(fileID, 'tStart = %gps\r\n', tStart/ps);
    fprintf(fileID, 'tFinish = %gps\r\n', tFinish/ps);
    fprintf(fileID, 'Propagation Method: %d\r\n', propagationMethod);

    fprintf(fileID, '\r\nUnits:\r\n');
    fprintf(fileID, '------\r\n');

    fprintf(fileID, 'hBar = %.16e\r\n', hBar);
    fprintf(fileID, 'c = %.16e\r\n', c);
    fprintf(fileID, 'eV = %.16e\r\n', eV);
    fprintf(fileID, 'amu = %.16e\r\n', amu);
    fprintf(fileID, 'A = %.16e\r\n', A);
    fprintf(fileID, 'ps = %.16e\r\n', ps);

    fprintf(fileID, '\r\nFull Accuracy:\r\n');
    fprintf(fileID, '--------------\r\n');
    fprintf(fileID, 'lx = %.16e \r\nly = %.16e \r\nlz = %.16e \r\n', lx, ly, lz);
    fprintf(fileID, 'dx = %.16e \r\ndy = %.16e \r\ndz = %.16e \r\n', dx, dy, dz);
    fprintf(fileID, 'dt = dt0 = %.16e\r\n', dt0);
    fprintf(fileID, 'tStart = %.16e\r\n', tStart);
    fprintf(fileID, 'tFinish = %.16e\r\n', tFinish);
    
    fprintf(fileID, '\r\n');
    
    fprintf(fileID, '\r\nPotential Setup\r\n');
    fprintf(fileID, '--------------\r\n');
    fprintf(fileID, 'alpha = %d\r\n', alphaIn);
    fprintf(fileID, 'xSigma = %g*A. ySigma = %g*A.\r\n', xSigmaIn/A, ySigmaIn/A);
    fprintf(fileID, 'gaussianPeakVal = %g*A\r\n', gaussPeakValIn/A);
    fprintf(fileID, 'wellDepth = %g*meV\r\n', wellDepthIn/(1e-3*eV));
    fprintf(fileID, 'saveLocation = %s\r\n', saveLocation);
    
    fprintf(fileID, '\r\n');
    
    fclose(fileID);
    
    % Save potential
    save('potentialUsed.mat', 'V');
end