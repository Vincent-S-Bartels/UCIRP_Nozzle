clear; clc; close all;
%% Nozzle Calculator Inputs
OFratio=input('Enter O/F wt. ratio (2.50:.05:3.50):');    %enter O/F wt ratio 2.50:.05:3.50
tRadius =input('Enter throat RADIUS in inches:');    % inches --later converted to metric
cPressure = input('[Chose 500 or 550 PSI]\nEnter Chamber Pressure in PSI: ');

%% Constants:
tRadius = tRadius*0.0254; % inches to meters
tArea = pi * tRadius^2; % m^2
%chosen because diameter of chamber should be ~3-4x throat crossection 
cArea = 3.5*tArea; % m^2
cRadius = sqrt(cArea/pi); % m
phi = 50; %pressure ratio

%% Calculated:
[massDot, gamma] = getMassFlowRate(OFratio, cPressure, tArea);
cPressure = cPressure * 6894.757; %psi to pascal
machExit = sqrt((phi^((gamma - 1)/gamma) - 1)*(2/(gamma - 1)) ); 
%uses isentropic flow relation for pressure and pressure ratio phi
epsilon = (1/machExit) * ((1 + ((gamma - 1)/2) * machExit^2)/...
    ((gamma+1)/2))^((gamma + 1)/(2*(gamma-1)));
%area ratio of exit to throat, based on isentropic flow relations
eArea = tArea * epsilon; % m^2
eRadius = sqrt(eArea/pi); % m
% eRadiusInches = eRadius*39.37008

%% Nozzle Design
% need a y(x), where y is the distance from centerline;
thetaInletCone = 15;
lInlet = (cRadius - tRadius)/tand(thetaInletCone); % meters
lThroat = (1/8)/39.37008; % meters
thetaOutletCone = 15; % Degrees
lOutlet = (2*tRadius * (sqrt(epsilon) - 1)) /(2*tand(thetaOutletCone)) ;
% maximum length based on 15 degree conical nozzle.
L1 = lInlet + lThroat;
LT = L1 + lOutlet;

%% Variables for crossection creation
x = linspace(0, (lInlet + lThroat + lOutlet), 2000);
y = zeros(length(x), 1)';

%for outlet portion, form of y = a + bx + cx^2 +dx^3
theta1 = -15;  % angle from chamber to inlet 
theta2 = 0;    % angle from inlet to throat, needs to be 0
theta3 = 5;    % angle from throat into outlet
theta4 = 15;   % angle from outlet to exit

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
j = 1;
A = [x', y'];
T = array2table(A, "VariableNames",{'Length (m)', 'Distance from Centerline (m)'});

writetable(T, 'NozzleDimensions.xlsx', 'sheet', j);

x = x*39.37008;
y = y*39.37008;
figure(1)
hold on
ylim([0,1.6])
plot(x,y)
title('Crossection of Nozzle')
xlabel('Length in Inches')
ylabel('Distance from Centerline in Inches')