function SetupBrownianMotionGaussians(displayAdsorbateAnimation, realTimePlotting)
    global A ps amu lx ly nx ny dx dy tStart tFinish gaussianPositions gaussianPositionsTimes numAdsorbates dt0 numIterations Browniefile savingBrownianPaths;
    xmax = [(lx - dx) (ly - dy)];
    vmax = 5*A/ps; % max initial speed in A/ps

    % Decay constant (for damped motion)
    gamma = 1*1e12;

    % Mass
    m = 3.0160293*amu;

    gaussianPositions = zeros(numAdsorbates, 2, numIterations+1);
    gaussianPositionsTimes = zeros(numAdsorbates, 2, numIterations+1);
    
    numDims = 2;
    if ny == 1
        numDims = 1;
    end
    for adsorbateNum = 1:numAdsorbates
        for dimNum = 1:numDims
            % Random initial starting conditions between 0 and xmax.
            x0 = xmax(dimNum)*rand;
            v0 = vmax*2*(rand - 0.5);

            y0 = [x0 v0];

            % Random motion parameters
            impulsesPerPs = 0.5;
            numImpulses = round(tFinish*impulsesPerPs/ps);

            % Random motion impulses at given times
            xi = 2*(rand(1, numImpulses) - 0.5);
            xi_t = linspace(tStart, tFinish+dt0, numImpulses);
      
            % Do integration
            [t, y] = ode45(@(t,y) RandomMotionDifferentialEquation(t,y, xi, xi_t, gamma, m), [tStart tFinish+dt0], y0);%adding dt0s due to overuse of interp1

<<<<<<< HEAD
            xtQuery = linspace(tStart, tFinish, numIterations+1);
=======
            xtQuery = linspace(tStart, tFinish, gaussianPositionsDiv);
>>>>>>> Added interpolation function for adsorbate position
            xQuery = interp1(t, y(:,1), xtQuery);
            
            gaussianPositions(adsorbateNum, dimNum, :) = xQuery;
            gaussianPositionsTimes(adsorbateNum, dimNum, :) = xtQuery;
        end
    end
    
    % Cyclic Boundary Conditions - make adsorbates loop back round to other side of simulation if they move past simulation boundary
    gaussianPositions(:,1,:) = mod(gaussianPositions(:,1,:), lx - dx);
    gaussianPositions(:,2,:) = mod(gaussianPositions(:,2,:), ly - dy);

    if(savingBrownianPaths)
       f=fopen(Browniefile, "wt");
       for i= 1:(1+numIterations)
           fprintf(f,"%f",(i-1)*dt0/ps);
           for ads = 1: numAdsorbates
               for dim = [1,2]
                   fprintf(f,"% f",gaussianPositions(ads,dim,i)/A);
               end
           end
           fprintf(f,"\n");
       end
       fclose(f);
    end

    figure %%%
    for adsorbateNum = 1:numAdsorbates
        plot(squeeze(gaussianPositions(adsorbateNum, 1, :))/A, squeeze(gaussianPositions(adsorbateNum, 2, :))/A, '.', 'MarkerSize', 1);
        hold on;
    end
    xlim([0 lx/A]);
    ylim([0 ly/A]);
    xlabel('lx /A');
    ylabel('ly /A');
    hold off
    
    % Display animated adsorbate motion before wavefunction propagation
    if displayAdsorbateAnimation
        figure
        DisplayAdsorbateAnimation(squeeze(gaussianPositions(:,1,:)), squeeze(gaussianPositions(:,2,:)), lx - dx, ly - dy);
    end
    
    % Prepare a new window if need to display real time animation
    if realTimePlotting
        figure
    end
end