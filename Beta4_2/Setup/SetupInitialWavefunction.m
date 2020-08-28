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
    vel = 800;             % m/s
    
    psi0 = zeros(size(V), 'gpuArray');
    
    psi0 = InitialiseGaussianWavefunction3D();
    %psi0 = InitialisePlaneWavefunction3D();
end