%% Documentation   
% 
%% Initialization
clear all; close all; clc; % addpath()
KT = karman_trefftz;


%% Parameters
% Main parameters affecting airfoil geometry
XC      = -0.2             ; % X-coord of cylinder center in Z-plane [m]
YC      = 0.1              ; % Y-coord of cylinder center in Z-plane [m]
tau_deg = 10               ; % Trailing edge [deg]
A       = 1                ; % Cylinder X-intersect in Z-plane [m] 
% Plotting parameters
n       = 100              ; % Number of points for airfoil

%% Derived parameters
l       = 2 - tau_deg / 180; % Karman-Trefftz "Lambda" parameter

% %% Airfoil shape
[xa, ya] = KT.shape(XC, YC, l, A, n);

%% Plot
figure()
plot(xa,ya)
