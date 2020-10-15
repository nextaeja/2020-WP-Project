% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function PrintEnergies()
    global mass vel;
    
    En = EnergyExpectation();
    Mom = MomentumExpectation();
    KE = sum(Mom.^2)/(2*mass);
    Vel = sqrt(2*KE/mass);
    percentDiffEnergy = 100*(En - KE)/En;
    diffVel = vel-Vel;
    percentDiffVel = 100*(vel - Vel)/vel;
    
    fprintf('En = %7d + %7di \nKE = %7d + %7di \nVel = %7d + %7di\nEnergy percent diff = %4d%% \nVel diff = %4d \nVel percent diff = %4d%% \n', real(En), imag(En), real(KE), imag(KE), real(Vel), imag(Vel), percentDiffEnergy, diffVel, percentDiffVel);
end