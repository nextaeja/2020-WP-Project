% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

randArray = 2*(rand(1, length(zOffset))-0.5);
error0 = 1e-3*max(zOffset);

error = error0*randArray;

zOffsetRand = error + zOffset;

%%
plot((x/A-100)*(6/5.5), zOffset/A/1.61, '-')
yyaxis right;
zError = abs((zOffset - zOffsetRand)./zOffset);
plot((x/A-100)*(6/5.5), zError, '-')
hold off
ylabel('Pointwise fractional error')

%%
avgGauss = mean(abs(zOffset));
zErrorConst = ((zOffset - zOffsetRand)/avgGauss);