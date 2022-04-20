clear all; clc; close all;
%% Inputs and constants
OFratio = linspace(2.5, 3.5, 21); % All OFratio in CEA
cPressure =  550;
cPressure = cPressure * 6894.757; %psi to pascal

tRadius =input('Enter throat RADIUS in inches:');    % inches --later converted to metric

tRadius = tRadius*0.0254; % inches to meters
tArea = pi * tRadius^2; % m^2
%chosen because diameter of chamber should be ~3-4x throat crossection
cArea = 3.5*tArea; % m^2
cRadius = sqrt(cArea/pi); % m
phi = 50; %pressure ratio

AngleforLength = [5, 15, 30, 45, 60]; % choosing to vary the enterance length
inletLength = (cRadius - tRadius)./tand(AngleforLength); % meters, triangle relation
inletThroat = (1/8)/39.37008; % meters
outletLength = (2*tRadius * (sqrt(epsilon) - 1)) /(2*tand(thetaOutletCone)); % meters, formula from M&ToP

%% Loop:
for i = 0: length(OFratio)
    [massDot, gamma] = getMassFlowRate(OFratio, cPressure, tArea);
    % returns the mass flow rate of the engine and average gamma between
    % throat and exit
    machExit = sqrt((phi^((gamma - 1)/gamma) - 1)*(2/(gamma - 1)) );
    %uses isentropic flow relation for pressure and pressure ratio phi
    epsilon = (1/machExit) * ((1 + ((gamma - 1)/2) * machExit^2)/...
        ((gamma+1)/2))^((gamma + 1)/(2*(gamma-1)));
    %area ratio of exit to throat, based on isentropic flow relations
    eArea = tArea * epsilon; % m^2
    eRadius = sqrt(eArea/pi); % m
    %% Nozzle Design

end