hBar = 1.054571800e-34; % Js
c = 299792458; % m/s
eV = 1.6021766208e-19; % J
amu =  1.660539040e-27; % kg
A = 1e-10; % m
ps = 1e-12; % s

% Setup
zCharacteristic = (1/2.06)*A;
sigmax = (5.50/6)*A;
zExtent = 1.61*A;
vel = 800;
mass = 3.0160293*amu;


for sigmaForward = [12:5:102]*A
    
    for lx = [100:50:2500]*A
        
        for lz = [100:50:2500]*A
            
            for nx = [512]
                
                for nz = [512]
                    % Nyquist
                    kMax = (mass*vel/hBar) + 5/sigmaForward;
                    nMin = max(kMax*lz/pi, kMax*lx/pi);

                    pxPerSigmax = (sigmax/lx)*nx;
                    pxPerZExtent = (zExtent/lz)*nz;
                    pxPerZCharacteristic = (zCharacteristic/lz)*nz;


                    % Criteria
                    % wavepacket must fit inside sim. 10sigmaForward <= lz
                    % px/zc etc criteria must be met
                    if(10*sigmaForward <= lz && pxPerSigmax >=1.6 && pxPerZExtent >=3.2 && pxPerZCharacteristic >= 0.77 && nMin < nx && nMin < nz)
                        fprintf('sigmaForward = %.2fA\t', sigmaForward/A);
                        fprintf('lx = %.2fA\t', lx/A);
                        fprintf('lz = %.2fA\t', lz/A);
                        fprintf('nMin = %.2f\t', nMin);
                        fprintf('nx = %.2f\t', nx);
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