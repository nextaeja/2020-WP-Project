% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

hBar = 1.054571800e-34; % Js
c = 299792458; % m/s
eV = 1.6021766208e-19; % J
amu =  1.660539040e-27; % kg
A = 1e-10; % m
ps = 1e-12; % s

% Setup
zCharacteristic = (1/2.06)*A;
sigmax = 2*(5.50/6)*A;
sigmay = 2*(5.50/6)*A;
zExtent = 2*1.61*A;
vel = 400;
mass = 3.0160293*amu;


for sigmaForward = [4:1:30]*A
    
    for lx = [20:10:300]*A
        
        for ly = [20:10:300]*A

            for lz = [20:10:300]*A

                for nx = [64 128 256 512]

                    for ny = [64 128 256 512]

                        for nz = [64 128 256 512 1024]
                            % Nyquist
                            kMax = (mass*vel/hBar) + 5/sigmaForward;
                            nMin = max([kMax*lx/pi, kMax*ly/pi, kMax*lz/pi]);

                            pxPerSigmax = (sigmax/lx)*nx;
                            pxPerSigmay = (sigmay/ly)*ny;
                            pxPerZExtent = (zExtent/lz)*nz;
                            pxPerZCharacteristic = (zCharacteristic/lz)*nz;


                            % Criteria
                            % wavepacket must fit inside sim. 10sigmaForward <= lz
                            % px/zc etc criteria must be met
                            if(10*sigmaForward <= lz &&...
                                    pxPerSigmax >=1.6 &&...
                                    pxPerSigmay >=1.6 &&...
                                    pxPerZExtent >=3.2 &&...
                                    pxPerZCharacteristic >= 2*0.77 &&...
                                    nMin < nx &&...
                                    nMin < ny &&...
                                    nMin < nz &&...
                                    nx == ny &&...
                                    nx*ny*nz < 256^3)
                                
                                fprintf('sigmaForward = %.2fA\t', sigmaForward/A);
                                fprintf('lx = %.2fA\t', lx/A);
                                fprintf('ly = %.2fA\t', ly/A);
                                fprintf('lz = %.2fA\t', lz/A);
                                fprintf('nMin = %.2f\t', nMin);
                                fprintf('nx = %.2f\t', nx);
                                fprintf('ny = %.2f\t', ny);
                                fprintf('nz = %.2f\t', nz);
                                fprintf('pxPerSigmaX = %.2f\t', pxPerSigmax);
                                fprintf('pxPerZExtent = %.2f\t', pxPerZExtent);
                                fprintf('pxPerZCharacteristic = %.2f\n', pxPerZCharacteristic);
                            end

                        end
                    end
                end
            end
        end
    end
end