% Input = 'decayType' = integer, values correspond to:
% 1 = Exponential
% 2 = Morse
% 3 = Morse-like, with alpha value specified by 'inParameterIfNeeded'
%
% Note: parameterIfNeeded is unused for decayType 1 and 2, so should be set arbitrarily to 0
%
function V = SetupDynamicGaussianPotential(sp, decayType, inParameterIfNeeded, xSigmaIn, ySigmaIn, gaussPeakValIn, wellDepthIn, x0, y0)
    global nx ny nz lx ly dx dy dz eV A numAdsorbates; % Needed in function
    % global V; % Set in function
    
    % Gaussian Properties
    %x0 = (lx - dx)/2;   % x-centre of Gaussian
    %y0 = (ly - dy)/2;   % y-centre of Gaussian
    xSigma = xSigmaIn;%2*(5.50/6)*SIUnits.A;       % x standard deviation
    ySigma = ySigmaIn;%2*(5.50/6)*SIUnits.A;       % y standard deviation
    gaussianPeakVal = gaussPeakValIn;%2*1.61*SIUnits.A;   % peak value of Gaussian
    
    % Setup 1D x, y, z grids
    xGrid1D = gpuArray(sp.dx*(0:sp.nx-1));
    yGrid1D = gpuArray(sp.dy*(0:sp.ny-1));
    zGrid1D = gpuArray(sp.dz*(0:sp.nz-1));
    
    % Setup 2D x and y grids for 2D Gaussian calculation
    xGrid2D = repmat(xGrid1D, sp.ny, 1);   % Created dimensions are [ny nx]
    xGrid2D = permute(xGrid2D, [2 1]);  % Permute to [nx ny]
    yGrid2D = repmat(yGrid1D, sp.nx, 1);   % Created [nx ny]. No permutation needed
    
    % Setup new blank 3D gpu array for V
    V = zeros(sp.nx, sp.ny, sp.nz, 'gpuArray');
    
    % Use loop to calculate for all adsorbates
    zOffset2D = zeros(sp.nx, sp.ny, 'gpuArray');
    xSigma2D(1:sp.nx, 1:sp.ny) = xSigma; %%THIS ISNT A GPU ARRAY
    ySigma2D(1:sp.nx, 1:sp.ny) = ySigma; %%THIS ISNT A GPU ARRAY
    for adsorbateNum = 1:sp.numAdsorbates
        % Setup 2D constants to use as input arguments for arrayfun use
        x02D(1:sp.nx, 1:sp.ny) = x0(adsorbateNum);
        y02D(1:sp.nx, 1:sp.ny) = y0(adsorbateNum); %%THESE AREN'T GPU ARRAYS

        % Calculate 2D then 3D Gaussian grid - i.e. zoffset grid
        zOffset2D = zOffset2D + gaussianPeakVal*arrayfun(@gaussian2DGrid, xGrid2D, yGrid2D, x02D, y02D, xSigma2D, ySigma2D);
    end
    zOffset3D = repmat(zOffset2D, 1, 1, sp.nz); % repeat matrix in (x, y, z) by (1, 1, nz)
    
    % Setup 3D z grid
    zGrid3D = repmat(zGrid1D, sp.nx, 1, sp.ny); % Created [nx nz ny]
    zGrid3D = permute(zGrid3D, [1 3 2]);  % Permuted to [nx ny nz]
    
    % Calculate zEffective 3D grid
    zEffective3D = zGrid3D - zOffset3D;
    
    % Create constants to pass to arrayfun GPU calculation of 1D z potential
    %Vmax = 100e-3*eV;
    zCharacteristic = (1/2.06)*SIUnits.A; % Vroot = pt. where V = 0 i.e. the z value where V = 0
    zOffset = 0*SIUnits.A; % Shift entire V away from boundary to stop Q.Tunneling through V
    
    wellDepth = wellDepthIn;%10e-3*eV;
    wellMinZPt = 2*SIUnits.A;
    
    % Set Vmax from zCharacteristic
    Vmax = wellDepth*exp(wellMinZPt/zCharacteristic);
    
    %Calculate V from each point in zEffective
    switch decayType
        case 1
            % Setup matrices to perform calculation with arrayfun GPU method
            Vmax3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = Vmax;
            Vroot3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = zCharacteristic;
            VOffset3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = zOffset;
            
            % Call function
            zEffective3D = arrayfun(@zExpPotential1D, zEffective3D, Vmax3D, Vroot3D, VOffset3D);
            
        case 2
            % Setup matrices to perform calculation with arrayfun GPU method
            a = (1/wellMinZPt)*log(1+sqrt(1+(Vmax/wellDepth))); % Makes V = Vmax at z = 0
            wellDepth3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = wellDepth;
            wellMinZPt3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = wellMinZPt;
            a3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = a;
            
            % Call function
            zEffective3D = arrayfun(@morsePotential, zEffective3D, wellDepth3D, wellMinZPt3D, a3D);
            
        case 3
            % Setup matrices to perform calculation with arrayfun GPU method
            alpha = inParameterIfNeeded;
            a = (1/wellMinZPt)*log((1/2)*(alpha+sqrt(alpha^2 + 4*(Vmax/wellDepth))));
            
            wellDepth3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = wellDepth;
            wellMinZPt3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = wellMinZPt;
            a3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = a;
            alpha3D(1:sp.nx, 1:sp.ny, 1:sp.nz) = alpha;
            
            % Call function
            zEffective3D = arrayfun(@morseLikePotential, zEffective3D, wellDepth3D, wellMinZPt3D, a3D, alpha3D);
    end
    
    V = zEffective3D;
    
end
function gaussianVal = gaussian2DGrid(x, y, x0, y0, xSigma, ySigma)
    % Gaussian exp arguments
    xArg = -(x-x0)^2/(2*xSigma^2);
    yArg = -(y-y0)^2/(2*ySigma^2);
    
    gaussianVal = exp(xArg + yArg);
end
function potentialVal = zLinearPotential1D(z, Vmax, zCharacteristic)
    % Implement linear potential that goes to 0 at some cutoff
    if(z < zCharacteristic)
        potentialVal = (-Vmax/zCharacteristic)*z + Vmax;
    else % z >= Vroot
        potentialVal = 0;
    end
end
function potentialVal = zExpPotential1D(z, Vmax, zCharacteristic, offset)
    potentialVal = Vmax*exp(-(1/zCharacteristic)*(z - offset));
end
function potentialVal = morsePotential(z, wellDepth, wellMinZPt, a)
    potentialVal = wellDepth*((1-exp(-a*(z-wellMinZPt)))^2 - 1);
end
function potentialVal = morseLikePotential(z, wellDepth, wellMinZPt, a, alpha)
    potentialVal = wellDepth*(exp(-2*1e-10*(z - wellMinZPt)) - alpha*exp(-a*(z - wellMinZPt)));
end
