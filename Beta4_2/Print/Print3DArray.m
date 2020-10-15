% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function Print3DArray(array, nx, ny, nz)
    for k=1:nz
        for i=1:nx
            for j=1:ny
                if isreal(array)
                    fprintf("%e ", array(i, j, k));
                else
                    fprintf("(%e + i*%e) ", real(array(i, j, k)), imag(array(i, j, k)));
                end
            end
            fprintf("\n");
        end
        fprintf("\n");
    end
end

