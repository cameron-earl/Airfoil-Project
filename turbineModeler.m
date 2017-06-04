% this models a rotors of a wind turbine
% it prints coords that can be exported to a spreadsheet
% spreadsheet data can be used to create 3d model
% this also creates 2.5d figure with depth indicated by color

clc;

% define known variables
NACA4 = 2408;           % 4 digit NACA code
totalRadius = 23.5;     % total rotor radius in meters
cLIFT = 1.5;            % coefficient of lift
lambda = 6;             % local tip speed ratio
bladeCount = 3;         % number of blades
chordModifier = .8;          % optimal chord ratio

n = 160; %count of 2d slices to print and plot
r = linspace(1,totalRadius,n); %distance from base of rotor

% define matrices as functions of r
lambda_r = (r / totalRadius) * lambda;
theta_r = (2/3) * atan(1 ./ lambda_r);
chord_r = ((8 * pi * r) / (bladeCount * cLIFT)) .* (1 - cos(theta_r));
chord_r = chord_r * chordModifier;

% use for loop and airfoil function to plot
figure
axis square
tit = strcat('\fontsize{16}Graph of NACA ', num2str(NACA4));
title(tit)
hold on
cc = hsv(numel(r)); % use pretty colors
for i=1:numel(r)
    disp('r = ')
    disp(i)
    [Xu,Yu,Xl,Yl] = airfoil(NACA4,chord_r(i),theta_r(i));
    plot (Xu,Yu,Xl,Yl,'color',cc(i,:))
end
hold off
