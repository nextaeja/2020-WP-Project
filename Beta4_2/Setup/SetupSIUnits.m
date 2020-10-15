% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function SetupSIUnits()
    global hBar c eV amu A ps ;
    
    hBar = 1.054571800e-34; % Js
    c = 299792458; % m/s
    eV = 1.6021766208e-19; % J
    amu =  1.660539040e-27; % kg
    A = 1e-10; % m
    ps = 1e-12; % s
    
end