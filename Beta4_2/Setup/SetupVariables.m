% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function CUDA_pointers = SetupVariables()
    global lx ly lz nx ny nz; % Needed to run function
    global dx dy dz kSquared numAdsorbates numIterations; % Set in function
    
    % Distance between each grid point. Units = m
    dx = lx/nx;
    dy = ly/ny;
    dz = lz/nz;
    
    % k-space setup
    kSquared = KSquared();
    
    CUDA_pointers = allocate_all_arrays(nx, ny, nz, numAdsorbates, numIterations);
end