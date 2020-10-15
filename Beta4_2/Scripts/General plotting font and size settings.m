% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

%% Plotting

xlabel('x/A')
ylabel('z/A')
legend('Gaussian Potential', 'Hann Potential')

set(findall(gcf,'-property','FontSize'),'FontSize',16);
set(findall(gca,'-property','FontSize'),'FontSize',16);

title({'\fontsize{24}Gaussian and Hann Potentials, peak height = 1.61A, centred at 100A','Gaussian sigma = (5.50/6)A, Hann w = 5.50A'});
%title({'\fontsize{24}Contour plot of psi in k-space','and ring along which form factor is plotted'});

set(findall(gcf,'-property','FontName'),'FontName', 'Serif');

%% Saving
savingName = '';
print(savingName, '-dpng', '-r600')