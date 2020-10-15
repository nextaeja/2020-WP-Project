% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

% NOTE: Specular peaks are removed BY DEFAULT.
%       Must comment-out specular section if spectular peak wanted

function plotAlongRing2DInterpolated(fnToPlot, ringRadius, angleStart, angleFinish, lx, ly, lz)
%=============SAMPLING================================%
    % Put fnToPlot on CPU so can do interpolation (MATLAB implementation doesn't work on GPU)
    fnToPlot = gather(fnToPlot);
    
    % Firstly, interpolate given fnToPlot to create finer gridpoints
    [nx, ny, nz] = size(fnToPlot);
    
    % Change in length per pixel in kx, ky, kz
    dkx = 2*pi/lx;
    dky = 2*pi/ly;
    dkz = 2*pi/lz;
    
    % k-space axes
    kxScale = dkx*[0:nx/2-1, 0, -nx/2+1:-1];
    kxScale = fftshift(kxScale);
    kxScale = kxScale(2:end-1);
    kyScale = dky*[1];%0:ny/2-1, 0, -ny/2+1:-1];
    %kyScale = fftshift(kyScale);
    %kyScale = kyScale(2:end);
    kzScale = dkz*[0:nz/2-1, 0, -nz/2+1:-1];
    kzScale = fftshift(kzScale);
    kzScale = kzScale(2:end-1);
    
    % 3D versions of k-space axes
    [KX, KY, KZ] = meshgrid(kxScale, kyScale, kzScale);
    
    % Setup sampling sizes
    
    nxNew = 20000;
    nyNew = 1;
    nzNew = 20000;
    
    dkxNew = dkx*(nx/nxNew);
    dkyNew = dky*(ny/nyNew);
    dkzNew = dkz*(nz/nzNew);
    
    kxScaleNew = dkxNew*[0:nxNew/2-1, 0, -nxNew/2+1:-1];
    kxScaleNew = fftshift(kxScaleNew);
    kxScaleNew = kxScaleNew(2:end-1);
    kyScaleNew = dkyNew*[1];%0:nyNew/2-1, 0, -nyNew/2+1:-1];
    %kyScaleNew = fftshift(kyScaleNew);
    %kyScaleNew = kyScaleNew(2:end);
    kzScaleNew = dkzNew*[0:nzNew/2-1, 0, -nzNew/2+1:-1];
    kzScaleNew = fftshift(kzScaleNew);
    kzScaleNew = kzScaleNew(2:end-1);
    
    % q = query, i.e. to be queried
    [KXq, KYq, KZq] = meshgrid(kxScaleNew, kyScaleNew, kzScaleNew);
    
    %fnToPlot = ifftshift(fnToPlot);
    fnToPlot = fnToPlot(2:end-1, 1, 2:end-1);
    
%=============REMOVE=SPECULAR=PEAK====================%
    % Must remove peak in full k-space before interpolation
    % since full 2D/3D k-space is interpolated
    
    % Get Specular peak centre
    xStart0 = -dkx*(length(fnToPlot(:,1))/2 + 1);
    zStart0 = -dkz*(length(fnToPlot(1,:))/2 + 1);
    xSpecPeak = ringRadius*cos(pi/2);
    zSpecPeak = ringRadius*sin(pi/2);
    nxSpecPeak = round((xSpecPeak - xStart0)/dkx);
    nzSpecPeak = round((zSpecPeak - zStart0)/dkz);
    nySpecPeak = 1;
    
    % Replace Spec.Peak (extended in nx and nz if wavefunction windowed in nx and nz)
    % with values taken from outside the peak (i.e. from an nx line a few pixels away from the peak)
    nxShift = 0; % 0 for plane wave as then kx Spec.Peak confined to 1px
    nzShift = 8; % Chosen arbitrarily to ensure all peak is sufficiently removed
    
    for nxi = nxSpecPeak - nxShift : nxSpecPeak + nxShift
        for nzi = nzSpecPeak - nzShift:nzSpecPeak + nzShift
            fnToPlot(nxi, nySpecPeak, nzi) = fnToPlot(nxSpecPeak - nxShift - 1, nySpecPeak, nzi);
        end
    end
    % Do Sampling
    fnToPlotSampled = interpn(squeeze(KX), squeeze(KZ), squeeze(fnToPlot), squeeze(KXq), squeeze(KZq), 'spline');
    
    fnToPlotSampled = squeeze(fnToPlotSampled);
    KXq = squeeze((KXq));
    
%=============GATHER=ANGLE=DATA=======================%
    % Get nx ny grid points for different angles
    
    % Angular fidelity = dTheta
    angleFidelity = (2*pi/360)*0.5; % 0.5 degree fidelity
    
    % x = length in 1 direction
    % y = length in 2 direction
    % r = ringRadius
    % theta = angle from +ve 1 direction anticlockwise
    
    % Initialise
    theta = angleStart:angleFidelity:angleFinish;
    nxList = zeros(length(theta), 1, 'gpuArray');
    nzList = zeros(length(theta), 1, 'gpuArray');
    fnValsToPlot = zeros(length(theta), 1, 'gpuArray');
    kxValsToPlot = zeros(length(theta), 1, 'gpuArray');
    
    xStart = -dkxNew*(length(fnToPlotSampled(:,1))/2 + 1);
    yStart = -dkzNew*(length(fnToPlotSampled(1,:))/2 + 1);
    
    % Get nx ny positions to plot in fnToPlot
    for i = 1:length(theta)
        x = ringRadius*cos(theta(i));
        y = ringRadius*sin(theta(i));
        nxList(i) = round((x - xStart)/dkxNew);
        nzList(i) = round((y - yStart)/dkzNew);
    end
    
    % Get values along ring in fnToPlot
    for i = 1:length(theta)
        %fprintf('nx = %d, ny = %d, fnToPlot = %f\n', nx(i), ny(i), fnToPlot(nx(i), ny(i)));
        fnValsToPlot(i) = abs(fnToPlotSampled(nxList(i), nzList(i)));
        kxValsToPlot(i) = KXq(nxList(i), nzList(i));
    end
    
    % Plot log10 k-space value and normalise peak to 1
    fnValsToPlot = log10(fnValsToPlot);
    maxValToPlot = max(fnValsToPlot);
    fnValsToPlot = fnValsToPlot/maxValToPlot;
    
    % Plot
    A = 1e-10; % m
    %plot(kxValsToPlot*A, fnValsToPlot, '-', 'LineWidth', 0.25, 'Color', [0, 0, 0]+0.5);
    plot(kxValsToPlot*A, fnValsToPlot, 'b-', 'LineWidth', 2);
end