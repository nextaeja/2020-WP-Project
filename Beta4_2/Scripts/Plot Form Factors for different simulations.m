% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

% Plotting form factors for different simulations

% Units
hBar = 1.054571800e-34; % Js
c = 299792458; % m/s
eV = 1.6021766208e-19; % J
amu =  1.660539040e-27; % kg
A = 1e-10; % m
ps = 1e-12; % s

% =============================NEED TO CHANGE DEPENDING ON SIMULATION SETUP:=======================================%

% Change directory to simulation folder to load psi_end data. (folder containing sim00001, etc...)
cd('E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\PhysicalSimulations\Gaussian Morse-like Quick z_c matched\VaryPotWellDepth')
addpath(genpath('E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\Toolbox\Testing\Beta3_2'))
    
folderContents = dir;

simStart = 1; % Sim to start plotting with
simEnd = 1; % Number of simulations
simDisplayList = [0]; % List of specific sims to display

mass = 3.0160293*amu;

vel = -800;

lx = 250*A;
ly = 250*A;
lz = 160*A;


% =============================NEED TO CHANGE DEPENDING ON SIMULATION SETUP:=======================================%

kExpectationVal = abs(mass*vel/hBar);

legendString = strings(simEnd,1);

%figure('units','normalized','outerposition',[0 0 1 1])
hold on;

j = 1;

for i = 1:length(folderContents)
    if(startsWith(folderContents(i).name, 'sim'))
        simName = folderContents(i).name;
        simNumString = extractAfter(simName, strfind(simName, 'sim')+2);
        simNum = str2double(simNumString);
        if( (simNum >= simStart && simNum <= simEnd) || ismember(simNum, simDisplayList) )
            cd(folderContents(i).name);

            load psi_end.mat; % loads into psi

            [nx, ny, nz] = size(psi);

            legendString(j) = strcat('nx = ', num2str(nx), '. nz = ', num2str(nz), '. sim', num2str(simNum));

            kx = (2*pi/lx)*[0:nx/2-1, 0, -nx/2+1:-1];
            kx = fftshift(kx);
            %kx = kx(2:end-1);
            psiInput = fftshift(fftn(psi));
            %psiInput = fftshift(psiInput);
            %psiInput = psiInput(2:end-1, 2:end-1);
            
            %kx = ifftshift(kx);
            %psiInput = ifftshift(psiInput);
            
            kx3D = repmat(kx, ny, 1, nz); % Created dimensions are [ny nx nz]
            kx3D = permute(kx3D, [2 1 3]); % Permute to get [nx ny nz]
            %kx3D = fftshift(kx3D);
            
            % 2D
            plotAlongRing2D(squeeze(psiInput), kExpectationVal, 0, pi, 2*pi/lx, 2*pi/lz, squeeze(kx3D))
            
            % 3D
            %plotAlongRing2D(squeeze(psiInput(:,64,:)), kExpectationVal, 0, pi, 2*pi/lx, 2*pi/lz, squeeze(kx3D))

            cd('..')
            j = j + 1;
        end
    end
end
j = j - 1; % j now = numSims

xlabel('kx /m^{-1}');
ylabel(['log10(Value along ring, |k| = ' num2str(kExpectationVal, '%e') ' m)']);
legend(legendString);

% set(findall(gcf,'-property','FontSize'),'FontSize',16);
% set(findall(gca,'-property','FontSize'),'FontSize',16);
% 
% title({'\fontsize{24}Form factors for varying nx and nz', 'Peak normalised to 1'});
% 
% set(findall(gcf,'-property','FontName'),'FontName', 'Serif');
%%
saveName = strcat('');
print(saveName, '-dpng', '-r600');


    