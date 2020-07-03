% NOTE: Specular peaks are removed BY DEFAULT.
%       Must comment-out specular section if spectular peak wanted

function plotAlongRing2D(fnToPlot, ringRadius, angleStart, angleFinish, dx, dy, plotAgainst)
    % Get nx ny grid points for different angles
    
    % Angular fidelity = dTheta
    angleFidelity = (2*pi/360)*0.5; % 0.5 degree fidelity
    
    % x = length in 1 direction
    % y = length in 2 direction
    % r = ringRadius
    % theta = angle from +ve 1 direction anticlockwise
    
    % Initialise
    theta = angleStart:angleFidelity:angleFinish;
    nx = zeros(length(theta), 1, 'gpuArray');
    ny = zeros(length(theta), 1, 'gpuArray');
    fnValsToPlot = zeros(length(theta), 1, 'gpuArray');
    kxValsToPlot = zeros(length(theta), 1, 'gpuArray');
    
    xStart = -dx*(length(fnToPlot(:,1))/2 + 1);
    yStart = -dy*(length(fnToPlot(1,:))/2 + 1);
    
    % Get nx ny positions to plot in fnToPlot
    for i = 1:length(theta)
        x = ringRadius*cos(theta(i));
        y = ringRadius*sin(theta(i));
        nx(i) = round((x - xStart)/dx);
        ny(i) = round((y - yStart)/dy);
    end
    
    % Get values along ring in fnToPlot
    for i = 1:length(theta)
        %fprintf('nx = %d, ny = %d, fnToPlot = %f\n', nx(i), ny(i), fnToPlot(nx(i), ny(i)));
        fnValsToPlot(i) = abs(fnToPlot(nx(i), ny(i)));
        kxValsToPlot(i) = plotAgainst(nx(i), ny(i));
    end
    
    % Remove Central Specular Peak
    % Located at centre of fnValsToPlot, therefore remove central point
    % fnValsToPlot should have odd # entries, therefore mid-point is length/2 rounded UP
    specPeakLocation = ceil(length(fnValsToPlot)/2);
    fnValsToPlot(specPeakLocation) = 0; % Makes log -> -Inf, not plotted in Matlab
    
    % Plot log10 k-space value and normalise peak to 1
    fnValsToPlot = log10(fnValsToPlot);
    maxValToPlot = max(fnValsToPlot);
    fnValsToPlot = fnValsToPlot/maxValToPlot;
    
    % Plot
    plot(kxValsToPlot, fnValsToPlot, '-');
end