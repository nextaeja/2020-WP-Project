%% Load data

% change directory to sim folder of simulation you want to plot

load psi_end.mat;

psift = squeeze(fftshift(fftn(psi)));

% to get nx and nz arrays, need to breakpoint at the end of "plotAlongRing2D" function then run from there

%% Contour plot

figure;
hold on;
normVal = 1/log10(max(max(max(abs(fnToPlotSampled)))));
contourf(rot90(log10(abs(squeeze(fnToPlotSampled)))*normVal, 3));
%contourf((log10(abs(squeeze(fnToPlot)))));

colorbar;
plot(nxList, nzList, 'r', 'LineWidth', 2)

%%
indices = find(log10(abs((fnToPlotSampled))) < -2.5);
fnToPlotSampled(indices) = 0;

%% Labelling figure
xlabel('\fontsize{16}nx');
ylabel('\fontsize{16}nz');
title({'\fontsize{24}Contour plot of psi in k-space','and ring along which form factor is plotted'});
%legend('\fontsize{16}log10 psi in k-space', '\fontsize{16}Form factor plotted around this ring');
set(gca, 'FontSize', 16);
pbaspect([1 1 1]);

%% Saving
saveName = 'Countour plot and form factor plotting ring - incorrect - gaussian repulsive sim1';
print(saveName, '-dpng', '-r600');
