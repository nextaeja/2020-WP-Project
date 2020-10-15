% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

close all
clear all

hBar = 1.054571800e-34; % Js
c = 299792458; % m/s
eV = 1.6021766208e-19; % J
amu =  1.660539040e-27; % kg
A = 1e-10; % m
ps = 1e-12; % s

alphaAll = [0:0.2:2];
alpha2 = [0, 2];

xSigma0 = (5.50/6)*A;       % x standard deviation
ySigma0 = (5.50/6)*A;       % y standard deviation

gaussPeakVal0 = 1.61*A;   % peak value of Gaussian

wellDepth0 = 10e-3*eV;
%{
 Vary adsorbate height - so gaussPeakVal
 Vary from 1* to 4* init value
 
 saveLocation = 'E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\PhysicalSimulations\Gaussian Morse-like Quick z_c matched\VaryAdsorbateHeight';
 for j = 1:length(alpha2)
     for i = 2.4:0.02:2.5
         
         cd('E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\Toolbox\Testing\Beta3_1b');
         WavepacketPropagation_beta3_1b(alpha2(j), xSigma0, ySigma0, i*gaussPeakVal0, wellDepth0, saveLocation);
         close all
         
     end
 end

 Vary adsorbate width - so xSigma AND ySigma
 Vary from 1* to 4*
 
 saveLocation = 'E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\PhysicalSimulations\Gaussian Morse-like Quick z_c matched\VaryAdsorbateWidth';
 for j = 1:length(alpha2)
     for i = 1:0.5:6
         
         cd('E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\Toolbox\Testing\Beta3_1b');
         WavepacketPropagation_beta3_1b(alpha2(j), i*xSigma0, i*ySigma0, gaussPeakVal0, wellDepth0, saveLocation);
         close all
         
     end
 end

 Vary adsorbate tot size - so xSigma, AND ySigma, AND gaussPeakVal
 Vary from 1* to 4*
 
 saveLocation = 'E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\PhysicalSimulations\Gaussian Morse-like Quick z_c matched\VaryAdsorbateTotSize';
 for i = 1:4
     for j = 1:length(alphaAll)
         
         cd('E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\Toolbox\Testing\Beta3_1b');
         WavepacketPropagation_beta3_1b(alphaAll(j), i*xSigma0, i*ySigma0, i*gaussPeakVal0, wellDepth0, saveLocation);
         close all
         
     end
 end

 Vary Potential Well Depth for setup with (2*h, 2*w)
 
 saveLocation = 'E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\PhysicalSimulations\Gaussian Morse-like Quick z_c matched\VaryPotWellDepth';
 for i = 1:15
     alpha0 = 2; % 2 = Morse
     j = 2; % 2*h, 2*w
     
     cd('E:\Shared OS folder\University\Work\Year 4\Project\Code\Matlab\Toolbox\Testing\Beta3_1b');
     WavepacketPropagation_beta3_1b(alpha0, j*xSigma0, j*ySigma0, j*gaussPeakVal0, i*1e-3*eV, saveLocation);
     close all
 end 
 %}
    




