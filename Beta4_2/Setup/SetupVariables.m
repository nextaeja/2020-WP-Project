function sp = SetupVariables(sp)
    global dx dy dz kSquared; % Set in function
    
    % Distance between each grid point. Units = m
    
    sp.dx = sp.lx/sp.nx;
    sp.dy = sp.ly/sp.ny;
    sp.dz = sp.lz/sp.nz;
    
    sp.kSquared = KSquared(sp);
    
    % Set up global variables to be compatible with old code
    dx=sp.dx;
    dy=sp.dy;
    dz=sp.dz;
    
    % k-space setup
    kSquared = sp.kSquared;
end