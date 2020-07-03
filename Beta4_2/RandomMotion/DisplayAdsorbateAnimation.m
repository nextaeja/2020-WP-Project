function DisplayAdsorbateAnimation(x, y, xmax, ymax)
    h = animatedline('MarkerFaceColor', 'b', 'Marker', '.', 'LineStyle', 'none');
    axis([0, xmax, 0, ymax]);
    
    showOnlyCurrentPosition = false;
    
    for i = 1:length(x)
        if showOnlyCurrentPosition
            clearpoints(h)
        end
        
        addpoints(h, mod(x(:, i), xmax), mod(y(:, i), ymax));
        drawnow
    end
end