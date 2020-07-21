function wp=SetupInitialWavefunction(V,sp)
    global sigmaForward sigmaPerp psix0 psiy0 psiz0 theta phi vel mass psi0; % Set in function
    
    % Wavepacket properties
    si = SIUnits;
    wp = WaveFunction;
    wp.psi0 = zeros(size(V), 'gpuArray');

    wp.psi0 = InitialiseGaussianWavefunction3D(sp, wp);
    
    % below are global variable implementations to be discarded
    sigmaForward = 9*si.A;     % Standard deviation of Gaussian envelope in the propagation direction. Units: m
    sigmaPerp = 40000*si.A;       % Std of gaussian envelope perpendicular to propagation direction. Units: m
    psix0 = 45*si.A;           % Wavefunction centre x start position
    psiy0 = 45*si.A;           % Wavefunction centre y start position
    psiz0 = 45*si.A;           % Wavefunction centre z start position
    theta = 0;              % Theta angle of wavepacket normal
    phi = 0;                % Phi angle of wavepacket normal
    mass = 3.0160293*si.amu;   % Units = kg
    vel = -800;             % m/s
    
    psi0=wp.psi0;
    %psi0 = InitialisePlaneWavefunction3D();
    %psi0 = InitialisePlaneWavefunction3D();
end