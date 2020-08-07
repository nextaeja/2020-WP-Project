%
% MOMENTUM EXPECTATION
%
% RETURNS p = [px, py, pz]
% px, py, pz are scalars
%
% <p> = integral(psi* p-op psi)
% where p-op is momentum operator = -i*hBar*grad
function pExpectation = MomentumExpectation()
    global psi hBar lx ly lz nx ny nz;
    
    kx = (2*pi/lx)*[0:nx/2-1, 0, -nx/2+1:-1];    
    ky = (2*pi/ly)*[0:ny/2-1, 0, -ny/2+1:-1];    
    kz = (2*pi/lz)*[0:nz/2-1, 0, -nz/2+1:-1];
    
    kx3D = repmat(kx, ny, 1, nz); % Created dimensions are [ny nx nz]
    kx3D = permute(kx3D, [2 1 3]); % Permute to get [nx ny nz]

    ky3D = repmat(ky, nx, 1, nz); % Created [nx ny nz]

    kz3D = repmat(kz, nx, 1, ny); % Created [nx nz ny]
    kz3D = permute(kz3D, [1 3 2]); % Permute to [nx ny nz]
    
    psiFT = fftn(psi);
    
    psiFTx = 1i*kx3D.*psiFT;
    pxPsi = -1i*hBar*ifftn(psiFTx);
    
    psiFTy = 1i*ky3D.*psiFT;
    pyPsi = -1i*hBar*ifftn(psiFTy);
    
    psiFTz = 1i*kz3D.*psiFT;
    pzPsi = -1i*hBar*ifftn(psiFTz);
    
    pxExpectation = sum(sum(sum(conj(psi).*pxPsi)));
    pyExpectation = sum(sum(sum(conj(psi).*pyPsi)));
    pzExpectation = sum(sum(sum(conj(psi).*pzPsi)));
    
    pExpectation = [pxExpectation, pyExpectation, pzExpectation];
end
