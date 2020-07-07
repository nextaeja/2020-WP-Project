function SaveSimulationRunning(time)
    global ps;
    saveName = strcat('figure_time', num2str(time/ps),".png");%%%
    print(saveName, '-dpng', '-r400');
end