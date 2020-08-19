function SetupSpecifiedPaths(displayAdsorbateAnimation, realTimePlotting)
       global pathfile numAdsorbates tStart tFinish ps A dx dy lx ly gaussianPositions gaussianPositionsTimes numIterations
       f=fopen(pathfile,'r');
       paths = fscanf(f,'%f',[1+2*numAdsorbates,Inf]);%you need to specify at least 2 points per adsorbate per dimension in the form t1 a1x a1y a2x a2y... t2 a1x... (i.e. time is its own entry) (units ps and ang)
       fclose(f);                                     %also the highest t must be >tFinish  
       t = paths(1,:)*ps;
       xtQuery = linspace(tStart, tFinish, numIterations+1);
       gaussianPositions = zeros(numAdsorbates, 2, numIterations+1);
       gaussianPositionsTimes = zeros(numAdsorbates, 2, numIterations+1);
       
       for absno=[1:numAdsorbates]
           for dimno=[0:1]
               gaussianPositionsTimes(absno, 1+dimno, :) = xtQuery;
               xQuery=interp1(t, paths(2*absno+dimno,:)*A, xtQuery);
               gaussianPositions(absno, dimno+1, :) = xQuery;
           end
       end
       figure
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