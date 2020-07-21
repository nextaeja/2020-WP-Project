classdef WaveFunction

    properties
        % Wavepacket properties
        sigmaForward = 9*SIUnits.A;     % Standard deviatoin of Gaussian envelope in the propagation direction. Units: m
        sigmaPerp = 40000*SIUnits.A;       % Std of gaussian envelope perpendicular to propagation direction. Units: m
        psix0 = 45*SIUnits.A;           % Wavefunction centre x start position
        psiy0 = 45*SIUnits.A;           % Wavefunction centre y start position
        psiz0 = 45*SIUnits.A;           % Wavefunction centre z start position
        theta = 0;              % Theta angle of wavepacket normal
        phi = 0;                % Phi angle of wavepacket normal
        mass = 3.0160293*SIUnits.amu;   % Units = kg
        vel = -800;             % m/s
        psi;

        psi0;
        time;

    end
end