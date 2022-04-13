function [mdot, gamma] = getMassFlowRate(OFratio,pressure_chamber, area_throat)

mach_throat = 1;    % mach
R = 8.3144598;  % Joules/(mol*Kelvin)

switch pressure_chamber
    case 500
        T = readtable('CEA_Proccessed/CEAParameters(500).xlsx');
    case 550
        T = readtable('CEA_Proccessed/CEAParameters(550).xlsx');
end
pressure_chamber = pressure_chamber * 6894.757;
switch OFratio
    case 2.50		
		temperature_chamber = T{1,14};  % Kelvin
		pressure_exit = T{1,13}*10^5;   % Pascal
		molarmass_chamber = T{1,8}; % g/mol
		molarmass_throat = T{1,9};  % g/mol
		molarmass_exit = T{1,10};   % g/mol	
        gamma_chamber = T{1,5};    % cp/cv
        gamma_throat = T{1,6}; % cp/cv
        gamma_exit = T{1,7};   % cp/cv
        rho = T{1,2}; %
    case 2.55		
		temperature_chamber = T{2,14};  % Kelvin
		pressure_exit = T{2,13}*10^5;   % Pascal
		molarmass_chamber = T{2,8}; % g/mol
		molarmass_throat = T{2,9};  % g/mol
		molarmass_exit = T{2,10};   % g/mol
        gamma_chamber = T{2,5};    % cp/cv
        gamma_throat = T{2,6}; % cp/cv
        gamma_exit = T{2,7};   % cp/cv
    case 2.60		
		temperature_chamber = T{3,14};  % Kelvin
		pressure_exit = T{3,13}*10^5;   % Pascal
		molarmass_chamber = T{3,8}; % g/mol
		molarmass_throat = T{3,9};  % g/mol
		molarmass_exit = T{3,10};   % g/mol
        gamma_chamber = T{3,5};    % cp/cv
        gamma_throat = T{3,6}; % cp/cv
        gamma_exit = T{3,7};   % cp/cv
    case 2.65		
		temperature_chamber = T{4,14};  % Kelvin
		pressure_exit = T{4,13}*10^5;   % Pascal
		molarmass_chamber = T{4,8}; % g/mol
		molarmass_throat = T{4,9};  % g/mol
		molarmass_exit = T{4,10};   % g/mol
        gamma_chamber = T{4,5};    % cp/cv
        gamma_throat = T{4,6}; % cp/cv
        gamma_exit = T{4,7};   % cp/cv
    case 2.70		
		temperature_chamber = T{5,14};  % Kelvin
		pressure_exit = T{5,13}*10^5;   % Pascal
		molarmass_chamber = T{5,8}; % g/mol
		molarmass_throat = T{5,9};  % g/mol
		molarmass_exit = T{5,10};   % g/mol
        gamma_chamber = T{5,5};    % cp/cv
        gamma_throat = T{5,6}; % cp/cv
        gamma_exit = T{5,7};   % cp/cv
    case 2.75		
		temperature_chamber = T{6,14};  % Kelvin
		pressure_exit = T{6,13}*10^5;   % Pascal
		molarmass_chamber = T{6,8}; % g/mol
		molarmass_throat = T{6,9};  % g/mol
		molarmass_exit = T{6,10};   % g/mol 
        gamma_chamber = T{6,5};    % cp/cv
        gamma_throat = T{6,6}; % cp/cv
        gamma_exit = T{6,7};   % cp/cv
    case 2.80		
		temperature_chamber = T{7,14};  % Kelvin
		pressure_exit = T{7,13}*10^5;   % Pascal
		molarmass_chamber = T{7,8}; % g/mol
		molarmass_throat = T{7,9};  % g/mol
		molarmass_exit = T{7,10};   % g/mol
        gamma_chamber = T{7,5};    % cp/cv
        gamma_throat = T{7,6}; % cp/cv
        gamma_exit = T{7,7};   % cp/cv
    case 2.85		
		temperature_chamber = T{8,14};  % Kelvin
		pressure_exit = T{8,13}*10^5;   % Pascal
		molarmass_chamber = T{8,8}; % g/mol
		molarmass_throat = T{8,9};  % g/mol
		molarmass_exit = T{8,10};   % g/mol 
        gamma_chamber = T{8,5};    % cp/cv
        gamma_throat = T{8,6}; % cp/cv
        gamma_exit = T{8,7};   % cp/cv
    case 2.90		
		temperature_chamber = T{9,14};  % Kelvin
		pressure_exit = T{9,13}*10^5;   % Pascal
		molarmass_chamber = T{9,8}; % g/mol
		molarmass_throat = T{9,9};  % g/mol
		molarmass_exit = T{9,10};   % g/mol 
        gamma_chamber = T{9,5};    % cp/cv
        gamma_throat = T{9,6}; % cp/cv
        gamma_exit = T{9,7};   % cp/cv
    case 2.95		
		temperature_chamber = T{10,14};  % Kelvin
		pressure_exit = T{10,13}*10^5;   % Pascal
		molarmass_chamber = T{10,8}; % g/mol
		molarmass_throat = T{10,9};  % g/mol
		molarmass_exit = T{10,10};   % g/mol
        gamma_chamber = T{10,5};    % cp/cv
        gamma_throat = T{10,6}; % cp/cv
        gamma_exit = T{10,7};   % cp/cv
    case 3.00		
		temperature_chamber = T{11,14};  % Kelvin
		pressure_exit = T{11,13}*10^5;   % Pascal
		molarmass_chamber = T{11,8}; % g/mol
		molarmass_throat = T{11,9};  % g/mol
		molarmass_exit = T{11,10};   % g/mol 
        gamma_chamber = T{11,5};    % cp/cv
        gamma_throat = T{11,6}; % cp/cv
        gamma_exit = T{11,7};   % cp/cv
     case 3.05		
		temperature_chamber = T{12,14};  % Kelvin
		pressure_exit = T{12,13}*10^5;   % Pascal
		molarmass_chamber = T{12,8}; % g/mol
		molarmass_throat = T{12,9};  % g/mol
		molarmass_exit = T{12,10};   % g/mol 
        gamma_chamber = T{12,5};    % cp/cv
        gamma_throat = T{12,6}; % cp/cv
        gamma_exit = T{12,7};   % cp/cv
    case 3.10		
		temperature_chamber = T{13,14};  % Kelvin
		pressure_exit = T{13,13}*10^5;   % Pascal
		molarmass_chamber = T{13,8}; % g/mol
		molarmass_throat = T{13,9};  % g/mol
		molarmass_exit = T{13,10};   % g/mol
        gamma_chamber = T{13,5};    % cp/cv
        gamma_throat = T{13,6}; % cp/cv
        gamma_exit = T{13,7};   % cp/cv
    case 3.15		
		temperature_chamber = T{14,14};  % Kelvin
		pressure_exit = T{14,13}*10^5;   % Pascal
		molarmass_chamber = T{14,8}; % g/mol
		molarmass_throat = T{14,9};  % g/mol
		molarmass_exit = T{14,10};   % g/mol
        gamma_chamber = T{14,5};    % cp/cv
        gamma_throat = T{14,6}; % cp/cv
        gamma_exit = T{14,7};   % cp/cv
    case 3.20		
		temperature_chamber = T{15,14};  % Kelvin
		pressure_exit = T{15,13}*10^5;   % Pascal
		molarmass_chamber = T{15,8}; % g/mol
		molarmass_throat = T{15,9};  % g/mol
		molarmass_exit = T{15,10};   % g/mol
        gamma_chamber = T{15,5};    % cp/cv
        gamma_throat = T{15,6}; % cp/cv
        gamma_exit = T{15,7};   % cp/cv
    case 3.25		
		temperature_chamber = T{16,14};  % Kelvin
		pressure_exit = T{16,13}*10^5;   % Pascal
		molarmass_chamber = T{16,8}; % g/mol
		molarmass_throat = T{16,9};  % g/mol
		molarmass_exit = T{16,10};   % g/mol
        gamma_chamber = T{16,5};    % cp/cv
        gamma_throat = T{16,6}; % cp/cv
        gamma_exit = T{16,7};   % cp/cv
    case 3.30		
		temperature_chamber = T{17,14};  % Kelvin
		pressure_exit = T{17,13}*10^5;   % Pascal
		molarmass_chamber = T{17,8}; % g/mol
		molarmass_throat = T{17,9};  % g/mol
		molarmass_exit = T{17,10};   % g/mol
        gamma_chamber = T{17,5};    % cp/cv
        gamma_throat = T{17,6}; % cp/cv
        gamma_exit = T{17,7};   % cp/cv
    case 3.35		
		temperature_chamber = T{18,14};  % Kelvin
		pressure_exit = T{18,13}*10^5;   % Pascal
		molarmass_chamber = T{18,8}; % g/mol
		molarmass_throat = T{18,9};  % g/mol
		molarmass_exit = T{18,10};   % g/mol
        gamma_chamber = T{18,5};    % cp/cv
        gamma_throat = T{18,6}; % cp/cv
        gamma_exit = T{18,7};   % cp/cv
    case 3.40		
		temperature_chamber = T{19,14};  % Kelvin
		pressure_exit = T{19,13}*10^5;   % Pascal
		molarmass_chamber = T{19,8}; % g/mol
		molarmass_throat = T{19,9};  % g/mol
		molarmass_exit = T{19,10};   % g/mol
        gamma_chamber = T{19,5};    % cp/cv
        gamma_throat = T{19,6}; % cp/cv
        gamma_exit = T{19,7};   % cp/cv
    case 3.45		
		temperature_chamber = T{20,14};  % Kelvin
		pressure_exit = T{20,13}*10^5;   % Pascal
		molarmass_chamber = T{20,8}; % g/mol
		molarmass_throat = T{20,9};  % g/mol
		molarmass_exit = T{20,10};   % g/mol
        gamma_chamber = T{20,5};    % cp/cv
        gamma_throat = T{20,6}; % cp/cv
        gamma_exit = T{20,7};   % cp/cv
    case 3.50		
		temperature_chamber = T{21,14}; % Kelvin
		pressure_exit = T{21,13}*10^5;  % Pascal
		molarmass_chamber = T{21,8};    % g/mol
		molarmass_throat = T{21,9}; % g/mol
		molarmass_exit = T{21,10};  % g/mol
        gamma_chamber = T{21,5};    % cp/cv
        gamma_throat = T{21,6}; % cp/cv
        gamma_exit = T{21,7};   % cp/cv
end			
%% Calculate the Gas Constants & Heat Capacity
R_chamber = (R/molarmass_chamber)*1000; %Joules/(Kg*K)
R_throat = (R/molarmass_throat)*1000; %Joules/(Kg*K)
R_exit = (R/molarmass_exit)*1000; %Joules/(Kg*K)

cp_chamber = (gamma_chamber*R_chamber)/(gamma_chamber-1); %Joules/(Kg*K)
cp_throat = (gamma_throat*R_throat)/(gamma_throat-1);%Joules/(Kg*K)
cp_exit = (gamma_exit*R_exit)/(gamma_exit-1);%Joules/(Kg*K)
cv_chamber = (R_chamber)/(gamma_chamber-1); %Joules/(Kg*K)
cv_throat = (R_throat)/(gamma_throat-1); %Joules/(Kg*K)
cv_exit = (R_exit)/(gamma_exit-1);
%% Calculate the Gas Constants & Heat Capacity
R_chamber = (R/molarmass_chamber)*1000; %Joules/(Kg*K)

%% Calculate Mass Flow Rate
mdot = ((area_throat*pressure_chamber*mach_throat)/...
    (sqrt(temperature_chamber*R_chamber)))*((sqrt(gamma_throat*...
    (1+((gamma_throat-1)/2)*mach_throat^2)^...
    (-(gamma_throat+1)/(gamma_throat-1))))); % mass flow rate through an orifice at various mach numbers (choked is mach_throat = 1)
gamma = (gamma_exit + gamma_throat )/2;
end