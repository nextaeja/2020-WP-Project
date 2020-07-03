function SetupVariables()
    global lx ly lz nx ny nz; % Needed to run function
    global dx dy dz kSquared; % Set in function
    
    % Distance between each grid point. Units = m
    dx = lx/nx;
    dy = ly/ny;
    dz = lz/nz;
    
    % k-space setup
    kSquared = KSquared();
end