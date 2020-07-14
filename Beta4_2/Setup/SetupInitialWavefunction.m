function SetupInitialWavefunction()
    global A amu V; % Needed in function
    global sigmaForward sigmaPerp psix0 psiy0 psiz0 theta phi vel mass psi0; % Set in function
    
    % Wavepacket properties
    sigmaForward = 9*si.A;     % Standard deviatoin of Gaussian envelope in the propagation direction. Units: m
    sigmaPerp = 40000*si.A;       % Std of gaussian envelope perpendicular to propagation direction. Units: m
    psix0 = 45*si.A;           % Wavefunction centre x start position
    psiy0 = 45*si.A;           % Wavefunction centre y start position
    psiz0 = 45*si.A;           % Wavefunction centre z start position
    theta = 0;              % Theta angle of wavepacket normal
    phi = 0;                % Phi angle of wavepacket normal
    mass = 3.0160293*si.amu;   % Units = kg
    vel = -800;             % m/s
    
    psi0 = zeros(size(V), 'gpuArray');
    
    psi0 = InitialiseGaussianWavefunction3D();
    
    wp = WaveFunction;
    
    psi0 = InitialiseGaussianWavefunction3D(sp);
    %psi0 = InitialisePlaneWavefunction3D();
    %psi0 = InitialisePlaneWavefunction3D();
end