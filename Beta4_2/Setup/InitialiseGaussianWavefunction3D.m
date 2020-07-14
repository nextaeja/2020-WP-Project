%
% INITIALISE GAUSSIAN ENVELOPE PLANE WAVE
%
% RETURNS psi, initialised as gaussian weighted plane wave.
% lpkt (packet length) is taken to = +-5sigma, i.e. lpkt = 10sigma.
% The gaussian is placed so that its centre is lpkt/2 away from the end of the confining box (lpkt/2 in from lz).
%
function psi = InitialiseGaussianWavefunction3D(sp)
    global nx ny nz sigmaForward sigmaPerp psix0 psiy0 psiz0; % Needed in function
    
    % For plane wave
    [x3D, y3D, z3D, kx3D, ky3D, kz3D] = Setupxyzkxkykz3D();
    
    % For gaussian - k0 = k (~momentum) normalised to 1. k1 and k2 are perpendicular to k0.
    [k0x3D, k0y3D, k0z3D, k1x3D, k1y3D, k1z3D, k2x3D, k2y3D, k2z3D] = Setupk0k1k2();
    
    psix03D(1:sp.nx, 1:ny, 1:nz) = psix0;
    psiy03D(1:nx, 1:ny, 1:nz) = psiy0;
    psiz03D(1:nx, 1:ny, 1:nz) = psiz0;
    
    sigmaForward3D(1:nx, 1:ny, 1:nz) = sigmaForward;
    sigmaPerp3D(1:nx, 1:ny, 1:nz) = sigmaPerp;
    
    %psiUncropped = arrayfun(@GaussianWavefn, x3D, y3D, z3D, kx3D, ky3D, kz3D, k0x3D, k0y3D, k0z3D, k1x3D, k1y3D, k1z3D, k2x3D, k2y3D, k2z3D, sigmaForward3D, sigmaPerp3D, psix03D, psiy03D, psiz03D);
    psiUncropped = arrayfun(@GaussianExpWavefn, x3D, y3D, z3D, kx3D, ky3D, kz3D, k0x3D, k0y3D, k0z3D, k1x3D, k1y3D, k1z3D, k2x3D, k2y3D, k2z3D, sigmaForward3D, sigmaPerp3D, psix03D, psiy03D, psiz03D);
    
    psi = psiUncropped;
    
	psi = psi / sqrt(sum(sum(sum(abs(psi).^2))));
end
function gaussianWavefn = GaussianWavefn(x, y, z, kx, ky, kz, k0x, k0y, k0z, k1x, k1y, k1z, k2x, k2y, k2z, sigmaForward, sigmaPerp, psix0, psiy0, psiz0)
    % Plane wave: uses x, y, z, kx, ky, kz
    gaussianWavefn = exp(1i*(kx*x + ky*y + kz*z));
    
    % Get positions relative to wavefn centre = psix0, psiy0, psiz0
    posDiffx = x - psix0;
    posDiffy = y - psiy0;
    posDiffz = z - psiz0;
    
    % Gaussian envelope in propagation direction (forward direction = k0 direction)
    gaussianWavefn = gaussianWavefn*exp(-0.5*((k0x*posDiffx + k0y*posDiffy + k0z*posDiffz)/sigmaForward)^2);
    
    % Gaussian envelope in perpendicular direction
    gaussianWavefn = gaussianWavefn*exp(-0.5*((k1x*posDiffx + k1y*posDiffy + k1z*posDiffz)/sigmaPerp)^2)*exp(-0.5*((k2x*posDiffx + k2y*posDiffy + k2z*posDiffz)/sigmaPerp)^2);
    
end
function gaussianExpWavefn = GaussianExpWavefn(x, y, z, kx, ky, kz, k0x, k0y, k0z, k1x, k1y, k1z, k2x, k2y, k2z, sigmaForward, sigmaPerp, psix0, psiy0, psiz0)
    % Plane wave: uses x, y, z, kx, ky, kz
    gaussianExpWavefn = exp(1i*(kx*x + ky*y + kz*z));
    
    % Get positions relative to wavefn centre = psix0, psiy0, psiz0
    posDiffx = x - psix0;
    posDiffy = y - psiy0;
    posDiffz = z - psiz0;
    
    % Gaussian envelope in propagation direction (forward direction = k0 direction)
    gaussianExpWavefn = gaussianExpWavefn*exp(-0.5*((k0x*posDiffx + k0y*posDiffy + k0z*posDiffz)/sigmaForward)^2);
    
    % MATRIX OPERATIONS NOT SUPPORTED INSIDE GPUARRAY FUNCTION
    % THEREFORE, DO ALL OPERATIONS COMPONENT-WISE BY HAND. VERY BUG PRONE.
    
    % Current position x,y,z = pos
    
    % Wavefunction centre = pos0 = psix0, psiy0, psiz0
    
    % Orthonormal basis vectors, k0 in wavefn propagation direction
    %k0 = [k0x k0y k0z]; % Not needed as Gaussian in k0 already implemented
    %k1 = [k1x k1y k1z];
    %k2 = [k2x k2y k2z];
    
    % Define shape of window, in units of sigmaPerp
    centralWindowHalf = 2*sigmaPerp;  % Size of central section = 1, from centre to outer edge
    falloff =           1*sigmaPerp;  % Size of outer section over which fn falls from 1 to 0
    
    % Work out posDiff components in k1 and k2 direction. Save in k1 and k2
    % k1 = (k1*posDiff)*k1;
    % k2 = (k2*posDiff)*k2;
    k1Component = k1x*posDiffx + k1y*posDiffy + k1z*posDiffz;
    k2Component = k2x*posDiffx + k2y*posDiffy + k2z*posDiffz;
%     
%     k1x = k1x*k1Component;
%     k1y = k1y*k1Component;
%     k1z = k1z*k1Component;
%     
%     k2x = k2x*k2Component;
%     k2y = k2y*k2Component;
%     k2z = k2z*k2Component;
%     
%     % Want a circle in perpendicular direction, therefore only consider k1 + k2 combination = kPerp
%     kPerpx = k1x + k2x;
%     kPerpy = k1y + k2y;
%     kPerpz = k1z + k2z;
%     
%     % Magnitude of current point away from wavefunction centre in direction perpendicular to propagation direction
%     kPerpMag = sqrt(kPerpx^2 + kPerpy^2 + kPerpz^2);
    
    kPerpMag = sqrt(k1Component^2 +k2Component^2);
    
%     if kPerpMag > sigmaPerp*10
%         gaussianExpWavefn = 1*gaussianExpWavefn;
%     else
%         gaussianExpWavefn = 0*gaussianExpWavefn;
%     end
    
    % Check to see if current point is in outermost region, where windowing fn = 0
    if(kPerpMag > centralWindowHalf + falloff)
        % Windowing fn = 0
        gaussianExpWavefn = 0*gaussianExpWavefn;
    
    % Else, check if current point is in central region, where windowing fn = 1
    elseif(kPerpMag < centralWindowHalf)
        % Windowing fn = 1
        gaussianExpWavefn = 1*gaussianExpWavefn;
        
    % Else, current point in falloff region of windowing fn, proportional to cos
    elseif(kPerpMag >= centralWindowHalf && kPerpMag <= centralWindowHalf + falloff)
        % Cos falloff dependence
        gaussianExpWavefn = gaussianExpWavefn*((cos(pi*(kPerpMag - centralWindowHalf)/sigmaPerp)+1)/2);
        
    % Else - something went wrong - print to user
    else
        % Print statements not allowed in GPUARRAY Functions...
        gaussianExpWavefn = -Inf - 1i*Inf;%NaN + 1i*NaN; % Force crash... Not ideal but better NaN than incorrect results
    end
    
    
end
