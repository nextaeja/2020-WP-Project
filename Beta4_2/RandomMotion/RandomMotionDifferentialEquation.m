% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function dydt = RandomMotionDifferentialEquation(t, y, xi, xi_t, gamma, m)

dydt = zeros(2,1);
xit = interp1(xi_t, xi, t);
dydt(1) = y(2);
dydt(2) = -gamma*y(2) + 500*gamma*xit;
%fprintf('%eps\n', t/1e-12);