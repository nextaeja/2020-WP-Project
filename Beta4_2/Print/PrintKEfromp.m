% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

% KE from expectation value of p. Take this, and mod square to get p^2
function KE = PrintKEfromp()
    global mass;
    
    Mom = MomentumExpectation();
    KE = sum(Mom.^2)/(2*mass);
    
    fprintf('KE from <p>^2 = %.16e', KE);
end