% k^2 matrix setup
% Note: k^2 NOT defined 2*pi/Lx... but 1/Lx...
% Used to get del^2 in Fourier space by doing k^2*fft(psi)
function kSquared = KSquared(sp)
    
    kx = (2*pi/sp.lx)*[0:sp.nx/2-1, 0, -sp.nx/2+1:-1];    
    ky = (2*pi/sp.ly)*[0:sp.ny/2-1, 0, -sp.ny/2+1:-1];    
    kz = (2*pi/sp.lz)*[0:sp.nz/2-1, 0, -sp.nz/2+1:-1];
    

    % Most efficient method would use "repmat" to create nx*ny*nz matrices of repeating kx,ky,kz vectors
    % Then do kx3D.*kx3D and sum to get kSquared

    kx3D = repmat(kx, sp.ny, 1, sp.nz); % Created dimensions are [ny nx nz]
    kx3D = permute(kx3D, [2 1 3]); % Permute to get [nx ny nz]

    ky3D = repmat(ky, sp.nx, 1, sp.nz); % Created [nx ny nz]

    kz3D = repmat(kz, sp.nx, 1, sp.ny); % Created [nx nz ny]
    kz3D = permute(kz3D, [1 3 2]); % Permute to [nx ny nz]

    kSquared = kx3D.*kx3D + ky3D.*ky3D + kz3D.*kz3D;
    
    % Put in form of fft(psi) - i.e. -ve frequencies at the end of array
    %kSquared = ifftshift(kSquared); % - Already in fft setup from way kx,... initialised
    
end