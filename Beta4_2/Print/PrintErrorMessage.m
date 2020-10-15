% Copyright (c) 2020, Lorenzo Basso, Jack Lee, Matthew Zhang, Feiyang Chen
% Copyright (c) 2018, Francis Haghighi-Daly 
% All rights reserved.
% This file is part of the WooStOr - Wavepacket prOpopgatiOn using SpliT OperatR method, subject to the GNU/GPL-3.0-or-later.

function PrintErrorMessage(errMessage, notifyUserInTerminal)
    if(notifyUserInTerminal)
        fprintf('\nERROR\n');
        fprintf(errMessage);
        fprintf('\n\n');
    end
    
    fileID = fopen('InitilisationData.txt', 'a');
    fprintf(fileID, '\r\n===================\r\n');
    fprintf(fileID, 'ERROR\r\n');
    fprintf(fileID, '===================\r\n');
    fprintf(fileID, 'Error message:\r\n%s\r\n',errMessage);
    fclose(fileID);
end