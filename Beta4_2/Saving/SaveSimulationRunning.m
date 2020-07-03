function SaveSimulationRunning(time)
    global ps;
    saveName = strcat('figure_time', num2str(time/ps));
    print(saveName, '-dpng', '-r400');
end