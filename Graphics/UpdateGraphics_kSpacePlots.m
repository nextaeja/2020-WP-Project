function UpdateGraphics_kSpacePlots(currentTime, itNumCompleted) % itNum == Iteration Number
    
    global nx ny nz lz xScale zScale kxScale kzScale psi V eV A psi0
    
    % Adapted From Alderwick
    % 3 plots - always do 1D plot
    if nx ~= 1
        % Setup fft(psi) needed for k-space plots
        psiFT = fftshift(fftn(psi));
        
        % For k-space normalisation
        psi0FT = (fftn(psi0));
        
        %figure('units','normalized','outerposition',[0 0 1 1])
        
        if ny == 1
            % 1D:
            subplot(4,1,1); %2x2 grid, put this plot over positions 1 & 2
            xToPlot = round(nx/2);
            yToPlot = round(ny/2);
            zToPlot = round(nz/2);
    
            % Graphics scaling: psi normalised to 1.0 in full 3D space, so 1D amplitude small
            scaleFactorPsi = 1/max(abs(psi0(xToPlot, yToPlot, :)));
            scaleFactorPsiFT = 1/max(abs(psi0FT(xToPlot, yToPlot, :)));
            
            % Real-space plot
            % Plot wavefunction. Scale Axis on right - values 1 to -1
            yyaxis right;
            plot(zScale, scaleFactorPsi*squeeze(real(psi(xToPlot,yToPlot,:))), 'b-', zScale, scaleFactorPsi*squeeze(imag(psi(xToPlot,yToPlot,:))), 'r-'); %squeeze turns the 1*1*nz 3D-matrix to 1*nz 2D-matrix for plotting
            axis([0 lz -1 1]);
            xlabel('lz/m');
            ylabel('psi');
            title(['psi 1D, and V, along nx =' num2str(xToPlot) ', ny = ' num2str(yToPlot) ' against lz/m']);
            %Plot Potential. Scale axis on left.
            yyaxis left;
            plot(zScale,squeeze(V(xToPlot, yToPlot, :)));
            Vaxis = ylim; % Place origin through centre
            halfVRange = (Vaxis(2) - Vaxis(1))/2;
            axis([0 lz Vaxis(1)-halfVRange Vaxis(2)-halfVRange]);
            ylabel('V/J');
            
            % k-space plot
            % Plot wavefunction. Scale Axis on right - values 1 to -1
            subplot(4,1,2); %2x2 grid, put this plot over positions 1 & 2
            yyaxis right;
            shiftKzScale = fftshift(kzScale);
            shiftKxScale = fftshift(kxScale);
            plot(shiftKzScale(:), squeeze(abs(psiFT(xToPlot+1,yToPlot,:))), 'b-');%, shiftKzScale(:), squeeze(imag(psiFT(xToPlot+1,yToPlot,:))), 'r-'); %squeeze turns the 1*1*nz 3D-matrix to 1*nz 2D-matrix for plotting
            %axis([0 lz -1 1]);
            xlabel('kz/m');
            ylabel('psi');
            title(['k-space psi (psiFT) 1D, along nx =' num2str(xToPlot) ', ny = ' num2str(yToPlot) ' against kz/m^{-1}']);
            

            % 2D:
            subplot(4,1,[3,4]);
            surf(shiftKxScale(:), shiftKzScale(:), squeeze(abs(psiFT(:,yToPlot,:))));
            title('2D x-space plot of psi');
            shading interp;
            view(-20, 30);
            xlabel('z');
            ylabel('x');
            
            % Add time and step number to plot title.
            suptitle(['Time: ' num2str(currentTime, '%.16e') '. Step ' num2str(itNumCompleted) ' completed']);

        else
            % 1D:
            subplot(3,3,[1,2,3]);
            xToPlot = round(nx/2);
            yToPlot = round(ny/2);
            
            % Plot wavefunction. Scale Axis on right - values 1 to -1
            yyaxis right;
            plot(zScale, scaleFactorPsi*squeeze(real(psi(xToPlot, yToPlot,:))), 'b-', zScale, scaleFactorPsi*squeeze(imag(psi(xToPlot, yToPlot,:))), 'r-'); %squeeze turns the 1*1*nz 3D-matrix to 1*nz 2D-matrix for plotting
            axis([0 lz -1 1]);
            xlabel('lz/m');
            ylabel('psi');
            title(['psi 1D, and V, along nx =' num2str(xToPlot) ', ny = ' num2str(yToPlot) ' against lz/A']);
            
            %Plot Potential. Scale axis on left.
            yyaxis left;
            plot(zScale,squeeze(V(xToPlot, yToPlot,:)));
            Vaxis = ylim; % Place origin through centre
            halfVRange = (Vaxis(2) - Vaxis(1))/2;
            axis([0 lz Vaxis(1)-halfVRange Vaxis(2)-halfVRange]);
            ylabel('V/J');

            % 2D cut at 1/4 of ny:
            subplot(3,3,4);
            imagesc(zScale,xScale,squeeze(real(psi(:,floor(ny/4),:))));
            set(gca,'dataaspectratio',[1 1 1]);
            title('2D x-space slice at 1/4 of y');
            xlabel('z');
            ylabel('x');

            % 2D cut at 1/2 of ny:
            subplot(3,3,5);
            imagesc(zScale,xScale,squeeze(real(psi(:,floor(ny/2),:))));
            set(gca,'dataaspectratio',[1 1 1]);
            title('2D x-space slice at 1/2 of y');
            xlabel('z');
            ylabel('x');
            
            % 2D cut at 3/4 of ny:
            subplot(3,3,6);
            imagesc(zScale,xScale,squeeze(real(psi(:,floor(3*ny/4),:))));
            set(gca,'dataaspectratio',[1 1 1]);
            title('2D x-space slice at 3/4 of y');
            xlabel('z');
            ylabel('x');
            
%==============================================================================%
%======NEED=TO=CHECK=THESE=WORK=AND=KX=KY=AXES=LABELED=CORRECTLY===============%
%==============================================================================%
            
            % K-space 2D cut 1/4 of ny:
            subplot(3,3,7);
            imagesc(kzScale, kxScale,squeeze(real(psiFT(:,floor(ny/4),:))));
            set(gca,'dataaspectratio',[1 1 1]);
            title(['2D k-space slice at ny = ' num2str(floor(ny/4))]);
            xlabel('kz');
            ylabel('kx');
            
            % K-space 2D cut 1/2 of ny:
            subplot(3,3,8);
            imagesc(kzScale, kxScale,squeeze(real(psiFT(:,floor(ny/2),:))));
            set(gca,'dataaspectratio',[1 1 1]);
            title(['2D k-space slice at ny = ' num2str(floor(ny/2))]);
            xlabel('kz');
            ylabel('kx');
            
            % K-space 2D cut 3/4 of ny:
            subplot(3,3,9);
            imagesc(kzScale, kxScale,squeeze(real(psiFT(:,:,floor(3*nz/4)))));
            set(gca,'dataaspectratio',[1 1 1]);
            title(['2D k-space slice at nz = ' num2str(floor(3*nz/4))]);
            xlabel('kx');
            ylabel('ky');
            
            
            % Add time and step number to plot title. Use super title (suptitle)
            suptitle(['Time: ' num2str(currentTime, '%.16e') '. Step ' num2str(itNumCompleted) ' completed']);
            
        end
    end
    drawnow
end

