clear all; clc; close all;
%% Inputs and constants
OFratio=3.05;
%input('Enter O/F wt. ratio (2.50:.05:3.50):');
cPressure =  550;
cPressure = cPressure * 6894.757; %psi to pascal
tRadius =0.5; 
%input('Enter throat RADIUS in inches:');    % inches --later converted to metric

tRadius = tRadius*0.0254; % inches to meters
tArea = pi * tRadius^2; % m^2
%chosen because diameter of chamber should be ~3-4x throat crossection
cArea = 4*tArea; % m^2
cRadius = sqrt(cArea/pi); % m
phi = 50; %pressure ratio


[massDot, gamma] = getMassFlowRate(OFratio, 550, tArea);
% returns the mass flow rate of the engine and average gamma between
% throat and exit
machExit = sqrt((phi^((gamma - 1)/gamma) - 1)*(2/(gamma - 1)) );
%uses isentropic flow relation for pressure and pressure ratio phi
epsilon = (1/machExit) * ((1 + ((gamma - 1)/2) * machExit^2)/...
    ((gamma+1)/2))^((gamma + 1)/(2*(gamma-1)));
%area ratio of exit to throat, based on isentropic flow relations
eArea = tArea * epsilon; % m^2
eRadius = sqrt(eArea/pi); % m

AngleforLength = [35, 45, 60]; % choosing to vary the enterance length
inletLength = (cRadius - tRadius)./tand(AngleforLength); % meters, triangle relation
throatLength = (1/8)/39.37008; % meters
thetaOutletCone = 15;
outletLength = (2*tRadius * (sqrt(epsilon) - 1)) /(2*tand(thetaOutletCone)) * 0.8; % meters, formula from M&ToP
% 80% truncated cone nozzle.

theta1 = -15;  % angle from chamber to inlet
theta2 = 0;    % angle from inlet to throat, needs to be 0
theta3 = 0;    % angle from throat into outlet
theta4 = 0;   % angle from outlet to exit


%% Nozzle Design Loop:
for j = 1: length(AngleforLength)
    x = linspace(0, (inletLength(j) + throatLength + outletLength), 2000);
    y = zeros(length(x), 1)';
    L1 = inletLength(j) + throatLength;
    LT = L1 + outletLength;
    for i = 1:length(x)
        if x(i) < inletLength(j)
            m = tand(theta2) - tand(theta1);
            a = cRadius;
            b = tand(theta1);
            c = (3/inletLength(j)^2)*(tRadius - cRadius - tand(theta1)*inletLength(j) - (1/3)*inletLength(j)*m);
            d = m/(3*inletLength(j)^2) - 2*c/(3*inletLength(j));
            y(i) = a + b*(x(i)) + c*(x(i))^2 + d*(x(i))^3;
        elseif inletLength(j) <= x(i) && x(i) <(throatLength + inletLength(j))
            y(i) = tRadius;
        elseif (throatLength + inletLength(j)) <= x(i) && x(i) <= (inletLength(j) + throatLength + outletLength)
            m = tand(theta4) - tand(theta3);
            a = tRadius;
            b = tand(theta3);
            c = (3/outletLength^2)*(eRadius - tRadius - tand(theta3)*outletLength - (1/3)*outletLength*m);
            d = m/(3*outletLength^2) - 2*c/(3*outletLength);
            y(i) = a + b*(x(i) - L1) + c*(x(i) - L1)^2 + d*(x(i) - L1)^3;
        end

    end
    figure(j)
    hold on
    plot(x'*39.37008, y'*39.37008);
    title(AngleforLength(j));
    ylabel('Distance from Centerline');
    xlabel('Length of Nozzle');
    xlim([0,x(length(x))*39.37008]);
    ylim([0, 2])
    hold off
    A = [x', y'];
    T = array2table(A, "VariableNames",{'Length (m)', 'Distance from Centerline (m)'});
    thetaStr = string(AngleforLength(j));
    name = 'ASOP/InletAngle_' + thetaStr+'.xlsx';
    writetable(T, name);
end



