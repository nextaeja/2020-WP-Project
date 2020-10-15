% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function SetupWedgePotential()
    global nx ny nz eV A lz; % Needed in function
    global V; % Set in function

    % Set up wedge potential
    V = zeros(nx,ny,nz, 'gpuArray');
    Vmax = 15*eV; % Units = J
    
    % Define length of wedge
    %Vroot = floor(nz/5);    % Vroot = pt. where V = 0 i.e. the nz value where V = 0
    Vroot = nz*round(6*A/lz);
    % DIRECTLY FROM ALDERWICK:
    % Create linear potential: Have 50-1 array. Repeat for each pt. on x axis
    % Repeat this row for each pt. y axis
    Vrep = repmat(gpuArray([Vroot:-1:1]),nx,1, ny);
    Vrep = Vrep*Vmax;
    
    %Rearrange order because repmat puts them in a weird order
    Vrep = permute(Vrep, [1,3,2]);
    V(:, :, 1:Vroot) = Vrep;
end