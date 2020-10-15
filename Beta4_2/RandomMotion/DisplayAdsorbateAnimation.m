% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function DisplayAdsorbateAnimation(x, y, xmax, ymax)
    h = animatedline('MarkerFaceColor', 'b', 'Marker', '.', 'LineStyle', 'none');
    axis([0, xmax, 0, ymax]);
    
    showOnlyCurrentPosition = false;
    
    for i = 1:length(x)
        if showOnlyCurrentPosition
            clearpoints(h)
        end
        
        addpoints(h, mod(x(:, i), xmax), mod(y(:, i), ymax));
        drawnow
    end
end