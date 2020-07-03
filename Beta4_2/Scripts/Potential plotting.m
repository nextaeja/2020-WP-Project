%% Setup units
hBar = 1.054571800e-34; % Js
c = 299792458; % m/s
eV = 1.6021766208e-19; % J
meV = 1e-3*eV;
amu =  1.660539040e-27; % kg
A = 1e-10; % m
ps = 1e-12; % s

%% Hann Potential
z = [95:0.01:105]*A;
x = [95:0.01:105]*A;

x0 = 100*A;
z0 = 0.1*A%81.2*A;

h = 1.61*A;
w = 5.5*A;
Vt = 8*meV;
gamma = 2.06/A;

for i = 1:length(x)
    
    if(abs(x(i) - x0) < w/2)
        chi(i) = 0.5*h*(1+cos(2*pi*(x(i) - x0)/w));
    else
        chi(i) = 0;
    end
    
end

Vhann = zeros(length(x), length(z));

for j = 1:length(z)
    for i = 1:length(x)
        
        Vhann(i, j) = Vt*exp(-gamma*(z(j) - z0 - chi(i)));
        
    end
end

%% Gaussian potential
z = [95:0.1:105]*A;
x = [95:0.1:105]*A;

x0 = 100*A;                 % x-centre of Gaussian
xSigma = (5.50/6)*A;       % x standard deviation
gaussianPeakVal = 1.61*A;   % peak value of Gaussian

zOffset = gaussianPeakVal*exp(-0.5*((x-x0)/xSigma).^2);

Vmax = 100e-3*eV;
zCharacteristic = (1/2.06)*A;
gamma = 1/zCharacteristic;


Vgaussian = zeros(length(x), length(z));

for j = 1:length(z)
    for i = 1:length(x)
        
        Vgaussian(i, j) = Vmax*exp(-gamma*(z(j)-zOffset(i)));;
        
    end
end

%% Exponential potential

z = [0:0.1:10]*A;

Vmax = 100e-3*eV;
zCharacteristic = (1/2.06)*A; % Vroot = pt. where V = 0 i.e. the z value where V = 0
zOffset = 0*A; % Shift entire V away from boundary to stop Q.Tunneling through V
    
Vexp = Vmax*exp(-(1/zCharacteristic)*(z - zOffset));

%% Morse potential
x = [0:0.1:20]*A;
z = [0:0.1:20]*A;

Vmax = 10e-3*eV;
z0 = 5*A;
a1 = (1/z0)*log(1+sqrt(1+100e-3*eV/Vmax));

Vmorse = Vmax*((1-exp(-a1*(z-z0))).^2 - 1);

%% Morse-like morphing potential
z = [0:0.05:10]*A;

De = 10e-3*eV;
z0 = 2*A;
V0 = 100e-3*eV;

alpha = [0:0.2:2];

a = (1/z0)*log(0.5*(alpha+sqrt(alpha.^2 + 4*(V0/De))));

legendString = strings(length(a),1);

V = zeros(length(a), length(z));

for i = 1:length(alpha)
    V(i, :) = De*(exp(-2*a(i)*(z-z0)) - alpha(i)*exp(-a(i)*(z-z0)));
    
    legendString(i) = strcat('\alpha = ', num2str(alpha(i), '%.1f'));
end

% Plot all potentials for different alpha values
hold on
for i = 1:length(alpha);
    plot(z/A, V(i,:)/meV);
    xlabel('z/A')
    ylabel('V/meV')
    legend(legendString)
end

%% Plotting 1 - setup figure correctly before plotting so axes ratios end up correct
% If you change fontsize after plotting, doesn't replace numbers that don't fit, just shifts them up...

set(findall(gcf,'-property','FontSize'),'FontSize',24);
set(findall(gca,'-property','FontSize'),'FontSize',24);

%pbaspect([2 1 1])
set(gcf,'units','points','position',[600,100,600,600])

set(findall(gcf,'-property','FontName'),'FontName', 'Serif');

hold on

%% Plotting 2
%xlabel('log_{10}(# Steps)')
xlabel('k$_{z}$/\AA$^{-1}$','Interpreter','latex')
ylabel('log_{10}[mean |normalised error|]','Interpreter','tex')

%legend('\alpha = 0.0', '\alpha = 2  .0')
%legend(alphaLegendString, 'FontSize', 16)
%legend('Reciprocal-space \psi', 'Form Factor Plotting line')
legend('2D', '3D')

set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'));
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
grid minor

%title('log_{10}[Form Factors] - Normalised')
%title({'\fontsize{24}Morse Potential changing into Repulsive Exponential', 'by multiplying 2nd exponential by a factor \alpha and adjusting remaining parameters to keep constant potential height','Setup: Well Depth = 10meV, Well min at z = 2A. \alpha = 2 is Morse Potential'});
%title({'\fontsize{24}Contour plot of psi in k-space','and ring along which form factor is plotted'});


%% Saving
savingName = '- for report';
print(savingName, '-dpng', '-r600')