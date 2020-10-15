% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function SaveSimulationRunning(time)
    global ps;
    saveName = strcat('figure_time', num2str(time/ps),".png");%%%
    print(saveName, '-dpng', '-r400');
end