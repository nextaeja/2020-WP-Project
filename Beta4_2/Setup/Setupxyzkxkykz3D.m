% Creates x3D, ..., kx3D, ... matrices for other methods
% Used instead of for loops for quicker calculations
function [x3D, y3D, z3D, kx3D, ky3D, kz3D] = Setupxyzkxkykz3D(sp,wp)
    
    % Most efficient initialisatin methods use "repmat" to create sp.nx*sp.ny*sp.nz matrices of repeating x,y,z vectors
    % Put all matrices in gpu
    x3D = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    y3D = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    z3D = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    kx3D = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    ky3D = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    kz3D = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    
    x = sp.dx*(0:sp.nx-1);
    y = sp.dy*(0:sp.ny-1);
    z = sp.dz*(0:sp.nz-1);
    
    x3D = repmat(x, sp.ny, 1, sp.nz); % Created dimensions are [sp.ny sp.nx sp.nz]
    x3D = permute(x3D, [2 1 3]); % Permute to get [sp.nx sp.ny sp.nz]

    y3D = repmat(y, sp.nx, 1, sp.nz); % Created [sp.nx sp.ny sp.nz]

    z3D = repmat(z, sp.nx, 1, sp.ny); % Created [sp.nx sp.nz sp.ny]
    z3D = permute(z3D, [1 3 2]); % Permute to [sp.nx sp.ny sp.nz]
    
    k = wp.mass*wp.vel/SIUnits.hBar;
    kx = k*sin(wp.theta)*cos(wp.phi);
    kx3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = kx;

    ky = k*sin(wp.theta)*sin(wp.phi);
    ky3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = ky;

    kz = k*cos(wp.theta);
    kz3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = kz;
end