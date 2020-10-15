% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function SetupGraphicsVariables()
    global nx ny nz lx ly lz dx dy dz gfxSteps; % Needed in function
    global xScale yScale zScale kxScale kyScale kzScale; % Set in function
    
    % Only initialise xScale, yscale, zScale if used
    if gfxSteps > 0
        xScale = gpuArray(dx*(0:nx-1));
        yScale = gpuArray(dy*(0:ny-1));
        zScale = gpuArray(dz*(0:nz-1));
        
        kxScale = (2*pi/lx)*[0:nx/2-1, 0, -nx/2+1:-1];
        kyScale = (2*pi/ly)*[0:ny/2-1, 0, -ny/2+1:-1];    
        kzScale = (2*pi/lz)*[0:nz/2-1, 0, -nz/2+1:-1];

    else
        % Not worth effort creating arrays if not used
        xScale = 0;
        yScale = 0;
        zScale = 0;
        
        kxScale = 0;
        kyScale = 0;
        kzScale = 0;
    end
end