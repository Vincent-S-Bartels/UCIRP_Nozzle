clear; clc; close all;
%% Nozzle Calculator Inputs
OFratio=input('Enter O/F wt. ratio (2.50:.05:3.50):');    %enter O/F wt ratio 2.50:.05:3.50
tRadius =input('Enter throat RADIUS in inches:');    % inches --later converted to metric
cPressure = input('[Chose 500 or 550 PSI]\nEnter Chamber Pressure in PSI: ');

%% Constants:
tRadius = tRadius*0.0254;
tArea = pi * tRadius^2;
%chosen because diameter of chamber should be ~3-4x throat radius
cArea = 3.5*tArea;
cRadius = sqrt(cArea/pi);
phi = 50; %pressure ratio

%% Calculated:
[massDot, gamma] = getMassFlowRate(OFratio, cPressure, tArea);
cPressure = cPressure * 6894.757; %psi to pascal
machExit = sqrt((phi^((gamma - 1)/gamma) - 1)*(2/(gamma - 1)) );
epsilon = (1/machExit) * ((1 + ((gamma - 1)/2) * machExit^2)/...
    ((gamma+1)/2))^((gamma + 1)/(2*(gamma-1)));
eArea = tArea * epsilon;
eRadius = sqrt(eArea/pi);
eRadiusInches = eRadius*39.37008
%% Nozzle Design
% need a y(x), where y is the distance from centerline;
theta1 = -30;
theta2 = 0; % Degrees


thetaInletCone = 15;
lInlet = (cRadius - tRadius)/tand(thetaInletCone); % meters
lThroat = (1/8)/39.37008; % meters
thetaOutletCone = 15; % Degrees
lOutlet = (2*tRadius * (sqrt(epsilon) - 1)) /(2*tand(thetaOutletCone)) ;
% maximum length based on 15 degree conical nozzle.

%% Variables for crossection creation
x = linspace(0, (lInlet + lThroat + lOutlet), 2000);
y = zeros(length(x), 1)';

%for outlet portion, form of y = a + bx + cx^2 +dx^3

theta3 = 0; %% throat to outlet
theta4 = 0; %% outlet
L1 = lInlet + lThroat;
LT = L1 + lOutlet;

for i = 1:length(x)
    if x(i) < lInlet
        m = tand(theta2) - tand(theta1);
        a = cRadius;
        b = tand(theta1);
        c = (3/lInlet^2)*(tRadius - cRadius - tand(theta1)*lInlet - (1/3)*lInlet*m);
        d = m/(3*lInlet^2) - 2*c/(3*lInlet);
        y(i) = a + b*(x(i)) + c*(x(i))^2 + d*(x(i))^3;
    elseif lInlet <= x(i) && x(i) <(lThroat + lInlet)
        y(i) = tRadius;
    elseif (lThroat + lInlet) <= x(i) && x(i) <= (lInlet + lThroat + lOutlet)
        m = tand(theta4) - tand(theta3);
        a = tRadius;
        b = tand(theta3);
        c = (3/lOutlet^2)*(eRadius - tRadius - tand(theta3)*lOutlet - (1/3)*lOutlet*m);
        d = m/(3*lOutlet^2) - 2*c/(3*lOutlet);
        y(i) = a + b*(x(i) - L1) + c*(x(i) - L1)^2 + d*(x(i) - L1)^3;
    end
end
x = x*39.37008;
y = y*39.37008;
figure(1)
hold on
ylim([0,1.6])
plot(x,y)
title('Crossection of Nozzle')
xlabel('Length in Inches')
ylabel('Distance from Centerline in Inches')