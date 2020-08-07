function SaveSimulationEndData(simTimeEnd, simIterationEnd)
    global startingDirectory dt psi ps; % Needed in function
    
    % Append InitilisationData.txt file with "Simulation Finished". 'a' argument means append data (don't wipe).
    fileID = fopen('InitilisationData.txt', 'a');
    fprintf(fileID, '\r\n===================\r\n');
    fprintf(fileID, 'Simulation Finished\r\n');
    fprintf(fileID, '===================\r\n\r\n');
    fprintf(fileID, 'tFinal = %.16e\r\n', simTimeEnd);
    fprintf(fileID, 'dt used = %.16e\r\n', dt);
    fprintf(fileID, 'Number of steps = %d\r\n', simIterationEnd);
    fprintf(fileID, 'Compute time = %f', toc);
    fclose(fileID);

    % Save final psi to allow analysis later
    save('psi_end.mat', 'psi');
    
    %save('end_time.mat', 'simTimeEnd');
    
    saveName = strcat('figure_time', num2str(simTimeEnd/ps),".png");%%%
    print(saveName, '-dpng', '-r400');
    
    % Change back to starting directory - leave simlation saving location
    cd(startingDirectory);
end