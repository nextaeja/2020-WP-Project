% Creates x3D, ..., kx3D, ... matrices for other methods
% Used instead of for loops for quicker calculations
function [x3D, y3D, z3D, kx3D, ky3D, kz3D] = Setupxyzkxkykz3D()
    global mass vel theta phi nx ny nz dx dy dz hBar; % Needed in function
    
    % Most efficient initialisatin methods use "repmat" to create nx*ny*nz matrices of repeating x,y,z vectors
    % Put all matrices in gpu
    x3D = zeros(nx, ny, nz, 'gpuArray');
    y3D = zeros(nx, ny, nz, 'gpuArray');
    z3D = zeros(nx, ny, nz, 'gpuArray');
    kx3D = zeros(nx, ny, nz, 'gpuArray');
    ky3D = zeros(nx, ny, nz, 'gpuArray');
    kz3D = zeros(nx, ny, nz, 'gpuArray');
    
    x = dx*(0:nx-1);
    y = dy*(0:ny-1);
    z = dz*(0:nz-1);
    
    x3D = repmat(x, ny, 1, nz); % Created dimensions are [ny nx nz]
    x3D = permute(x3D, [2 1 3]); % Permute to get [nx ny nz]

    y3D = repmat(y, nx, 1, nz); % Created [nx ny nz]

    z3D = repmat(z, nx, 1, ny); % Created [nx nz ny]
    z3D = permute(z3D, [1 3 2]); % Permute to [nx ny nz]
    
    k = mass*vel/hBar;
    kx = k*sin(theta)*cos(phi);
    kx3D(1:nx, 1:ny, 1:nz) = kx;

    ky = k*sin(theta)*sin(phi);
    ky3D(1:nx, 1:ny, 1:nz) = ky;

    kz = k*cos(theta);
    kz3D(1:nx, 1:ny, 1:nz) = kz;
end