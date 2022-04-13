clear; clc; close all;
%% Nozzle Calculator Inputs
OFratio=input('Enter O/F wt. ratio (2.50:.05:3.50):');    %enter O/F wt ratio 2.50:.05:3.50
tRadius =input('Enter throat RADIUS in inches:');    % inches --later converted to metric
cPressure = input('[Chose 500 or 550 PSI]\nEnter Chamber Pressure in PSI: ');

%% Constants:
tRadius = tRadius*0.0254;
tArea = pi * tRadius^2;
cRadius = 3.5*tRadius; 
%chosen because diameter of chamber should be ~3-4x throat radius
cArea = cRadius^2 * pi;
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