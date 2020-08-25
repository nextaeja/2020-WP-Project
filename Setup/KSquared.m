% k^2 matrix setup
% Note: k^2 NOT defined 2*pi/Lx... but 1/Lx...
% Used to get del^2 in Fourier space by doing k^2*fft(psi)
function kSquared = KSquared()
    global nx ny nz lx ly lz; % Needed in function
    
    kx = (2*pi/lx)*[0:nx/2-1, 0, -nx/2+1:-1];    
    ky = (2*pi/ly)*[0:ny/2-1, 0, -ny/2+1:-1];    
    kz = (2*pi/lz)*[0:nz/2-1, 0, -nz/2+1:-1];
    

    % Most efficient method would use "repmat" to create nx*ny*nz matrices of repeating kx,ky,kz vectors
    % Then do kx3D.*kx3D and sum to get kSquared

    kx3D = repmat(kx, ny, 1, nz); % Created dimensions are [ny nx nz]
    kx3D = permute(kx3D, [2 1 3]); % Permute to get [nx ny nz]

    ky3D = repmat(ky, nx, 1, nz); % Created [nx ny nz]

    kz3D = repmat(kz, nx, 1, ny); % Created [nx nz ny]
    kz3D = permute(kz3D, [1 3 2]); % Permute to [nx ny nz]

    kSquared = kx3D.*kx3D + ky3D.*ky3D + kz3D.*kz3D;
    
    % Put in form of fft(psi) - i.e. -ve frequencies at the end of array
    %kSquared = ifftshift(kSquared); % - Already in fft setup from way kx,... initialised
    
end