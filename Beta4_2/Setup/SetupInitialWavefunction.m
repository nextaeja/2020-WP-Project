% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function SetupInitialWavefunction()
    global A amu V; % Needed in function
    global sigmaForward sigmaPerp psix0 psiy0 psiz0 theta phi vel mass psi0; % Set in function
    
    % Wavepacket properties
    sigmaForward = 9*A;     % Standard deviatoin of Gaussian envelope in the propagation direction. Units: m
    sigmaPerp = 40000*A;       % Std of gaussian envelope perpendicular to propagation direction. Units: m
    psix0 = 45*A;           % Wavefunction centre x start position
    psiy0 = 45*A;           % Wavefunction centre y start position
    psiz0 = 45*A;           % Wavefunction centre z start position
    theta = 0;              % Theta angle of wavepacket normal
    phi = 0;                % Phi angle of wavepacket normal
    mass = 3.0160293*amu;   % Units = kg
    vel = -800;             % m/s
    
    psi0 = zeros(size(V), 'gpuArray');
    
    psi0 = InitialiseGaussianWavefunction3D();
    %psi0 = InitialisePlaneWavefunction3D();
end