% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function CreateNewSimulationFolder()
    global savingDirectory; % Needed in function
    global serialNumber startingDirectory; % Set in function
    
    % Save starting directory to return to later
    startingDirectory = pwd;

    % Change directory to save simulation in designated folder
    if ~exist(savingDirectory, 'dir')
       mkdir(savingDirectory)
    end
    cd(savingDirectory);

    % Get folder contents to be able to create next simulation directory
    folderContents = dir;
    currentMaxSerial = 0;

    % Loop through all files in directory to find previous max. serial number
    for i = 1:length(folderContents)
        % Only look at simulation folders, named "sim..."
        if(startsWith(folderContents(i).name, 'sim'))
            currentSerialNum = extractAfter(folderContents(i).name, 'sim');
            currentSerialNum = str2double(currentSerialNum);
            if currentSerialNum > currentMaxSerial
                currentMaxSerial = currentSerialNum;
            end
        end
    end

    % Increment previous max. serial number to get current serial number
    serialNumber = currentMaxSerial + 1;

    % sprintf puts data into string format. %05d prints 5 digits with leading 0s if not long enough. If number longer than 5, prints all numbers.
    folderName = sprintf('sim%05d', serialNumber);

    mkdir(folderName); % Note: Need brackets otherwise will literally call folder "folderName".
    cd(folderName);
end