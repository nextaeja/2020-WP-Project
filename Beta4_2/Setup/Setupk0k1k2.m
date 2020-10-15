% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

% Creates normalised k vector (k0) and 2 perpendiculars (k1 and k2)
% Repeated in 3D matrices - used instead of for loops for quicker calculations
function [k0x3D, k0y3D, k0z3D, k1x3D, k1y3D, k1z3D, k2x3D, k2y3D, k2z3D] = Setupk0k1k2()
    global theta phi nx ny nz vel; % Needed in function
    
    % Most efficient initialisatin methods use "repmat" to create nx*ny*nz matrices of repeating x,y,z vectors
    % Put all matrices in gpu
    k0x3D = zeros(nx, ny, nz, 'gpuArray');
    k0y3D = zeros(nx, ny, nz, 'gpuArray');
    k0z3D = zeros(nx, ny, nz, 'gpuArray');
    k1x3D = zeros(nx, ny, nz, 'gpuArray');
    k1y3D = zeros(nx, ny, nz, 'gpuArray');
    k1z3D = zeros(nx, ny, nz, 'gpuArray');
    k2x3D = zeros(nx, ny, nz, 'gpuArray');
    k2y3D = zeros(nx, ny, nz, 'gpuArray');
    k2z3D = zeros(nx, ny, nz, 'gpuArray');
    
    % k vector - normalised to 1 automatically if using sin and cos functions as below
    kx = sin(theta)*cos(phi);
    ky = sin(theta)*sin(phi);
    kz = cos(theta);
    
    % Get perpendicular vectors
    k0 = sign(vel)*[kx ky kz]; % sign of velocity can change direction of k0
    kperp = null(k0(:).'); % Gets perpendicular vectors
    k1 = transpose(kperp(:,1)); % Transpose to get correct dims for later dot product
    k2 = transpose(kperp(:,2));
    
    % Fill 3D matrices
    k0x3D(1:nx, 1:ny, 1:nz) = k0(1);
    k0y3D(1:nx, 1:ny, 1:nz) = k0(2);
    k0z3D(1:nx, 1:ny, 1:nz) = k0(3);
    k1x3D(1:nx, 1:ny, 1:nz) = k1(1);
    k1y3D(1:nx, 1:ny, 1:nz) = k1(2);
    k1z3D(1:nx, 1:ny, 1:nz) = k1(3);
    k2x3D(1:nx, 1:ny, 1:nz) = k2(1);
    k2y3D(1:nx, 1:ny, 1:nz) = k2(2);
    k2z3D(1:nx, 1:ny, 1:nz) = k2(3);
end