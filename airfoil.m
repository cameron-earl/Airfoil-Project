function [Xu,Yu,Xl,Yl] = airfoil( NACA4, chordLength, theta)
%AIRFOIL generates airfoil shape with NACA4 designation and chord length
%   NACA4 is a 4 digit NACA designation for the airfoil. Input is a number.
%   chordLength is the total chord length. Input is a number.
%   theta is the number of radians it should be rotated counterclockwise
%   airfoil will be rotated around (.25chordLength,0)
%   i.e. airfoil(5512,1)

m = floor(NACA4 / 1000) / 100; % first digit is m, a percentage
% m = maximum camber in percentage of chord
p = floor(mod(NACA4,1000) / 100) / 10; % 2nd digit is p, a tenth value
%  p = position of the maximum camber in tenths of chord
t = mod(NACA4,100) / 100; % last two digits are t, a percentage
% t = maximum thickness of the airfoil in percentage of chord

numOfCoordinates = 10;
deltaX = 1 / numOfCoordinates; 

XA = [0:deltaX/10:p]; % equations are different before and after p
XB = [p:deltaX:1];
X = [XA, XB]; % concatenate arrays to full x domain for easy plotting

YcA = (m / p^2) * (2*p*XA - XA.^2);
YcB = (m / (1 - p)^2) * ((1 - 2*p) + 2*p*XB - XB.^2);
Yc = [YcA, YcB]; % y values for camber line

a = .2969;
b = -.126;
c = -.3516;
d = .2843;
e = -.1015;

YtA = (t / .2)*(a*XA.^.5 + b*XA + c*XA.^2 + d*XA.^3 + e*XA.^4);
YtB = (t / .2)*(a*XB.^.5 + b*XB + c*XB.^2 + d*XB.^3 + e*XB.^4);
% thickness distribution

thetaA = atan((2* m / p^2) * (p - XA));
thetaB = atan((2 * m / (1 - p)^2) * (p - XB));
% theta = arctan(dYc/dX)

XuA = XA - YtA .* sin(thetaA);
XuB = XB - YtB .* sin(thetaB);
Xu = [XuA, XuB]; % x values for upper edge of airfoil

YuA = YcA + YtA .* cos(thetaA);
YuB = YcB + YtB .* cos(thetaB);
Yu = [YuA, YuB]; % y values for upper edge of airfoil

XlA = XA + YtA .* sin(thetaA);
XlB = XB + YtB .* sin(thetaB);
Xl = [XlA, XlB]; % x values for lower edge of airfoil

YlA = YcA - YtA .* cos(thetaA);
YlB = YcB - YtB .* cos(thetaB);
Yl = [YlA, YlB]; % y values for lower edge of airfoil

% Plot unrotated, unsized airfoil
%figure
%axis equal
%hold on
%plot(Xu,Yu)
%plot(Xl,Yl)

% shape fully calculated, now multiplying all array elements by chordLength
%   to achieve proper size
X  = X  .* chordLength;
Yc = Yc .* chordLength;
Xu = Xu .* chordLength;
Yu = Yu .* chordLength;
Xl = Xl .* chordLength;
Yl = Yl .* chordLength;

% Now to rotate theta° around (.25chordLength,0)
% first translate x coordinates left
Xu = Xu - (.25 * chordLength);
Xl = Xl - (.25 * chordLength);

% create xy matrix for each line to aid in calculations
upper = [Xu;Yu];
lower = [Xl;Yl];

% define rotational matrix
R = [cos(theta),-sin(theta);sin(theta),cos(theta)];

% perform rotation
upper = R * upper;
lower = R * lower;

% extract x and y again
Xu = upper(1,:);
Yu = upper(2,:);
Xl = lower(1,:);
Yl = lower(2,:);

% translate x coordinates back
Xu = Xu + (.25 * chordLength);
Xl = Xl + (.25 * chordLength);

% plot the rotated airfoil
%plot(Xu,Yu)
%plot(Xl,Yl)
%hold off

% display a list of coordinates for excel, rotate to print easier

Xu = rot90(Xu)
Yu = rot90(Yu)
Xl = rot90(Xl)
Yl = rot90(Yl)
% return matrices to proper dimensions
Xu = rot90(Xu);
Yu = rot90(Yu);
Xl = rot90(Xl);
Yl = rot90(Yl);
end