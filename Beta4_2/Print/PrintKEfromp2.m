% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

% KE from expectation value of p^2. i.e. do *-k^2 rather than *1i*k in expectation p.
function KE = PrintKEfromp2()
    En = EnergyExpectation();
    fprintf('KE from <p^2> = %.16e', En);
    KE = En;
end