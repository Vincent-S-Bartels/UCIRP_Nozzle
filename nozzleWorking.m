clear; clc; close all;
%% Nozzle Calculator Knowns
tRadius = 0.5*0.0254; %% May change, check with Brian (inches to meters)
tArea = pi * tRadius^2;
cRadius = 3.5*tRadius; 
%chosen because diameter of chamber should be ~3-4x throat radius

disp('Hello World')
