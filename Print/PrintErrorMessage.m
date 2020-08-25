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