InitClear
% setFigurePath('/work/publications/phdthesis/figs_airfoil/');
% setMatFigPath('/work/publications/phdthesis/figs_airfoil/matfig/');
setFigureTitle(0);
setFigurePath('./');
setMatFigure(0);
setFigureFont('20');



%% Initialization
% require('PROFILES','v01')
addpath('../')


%% Parameters
% inputs for fully separated polar (in degrees)
alpha_merge=35;
delta_alpha=10;
dclda_fs=pi*(pi/180); % deg



%% 

% Loading profile data
% Data=load(['../data/tjaere11_ds.dat']);
% alpha=Data(:,1);
% Cl=Data(:,2);
% Cd=Data(:,3);
% Cl_inv=Data(:,5);
% Cl_fs=Data(:,6);
% f_st=Data(:,4);

Data=load(['data/FFA-W3-241-Re12e6.dat']);
alpha=Data(:,1);
Cl=Data(:,2);
Cd=Data(:,3);
Cm=Data(:,4);
%% Polar only
figure, hold all;box on; grid on;
plot(alpha,Cl,'k')
% [Cl_inv Cl_inv_sin alpha0]      = fPolarInviscid(alpha,Cl)                                        ; 
% plot(alpha,Cl_inv,'k--')
% vline(alpha0,'Color','k','LineStyle','-')

xlabel('\alpha [^o]')
ylabel('C_l [-]')
title('AirfoilFFAW3Cl');
xlim([-30 30])
ylim([-1.5 2])
%%


figure,hold all; box on ; grid on;
plot(alpha,Cd,'k')
xlabel('\alpha [^o]')
ylabel('C_d [-]')
title('AirfoilFFAW3Cd');
xlim([-30 30])
figure,hold all, box on; grid on;
plot(alpha,Cm,'k')
xlabel('\alpha [^o]')
ylabel('C_m [-]')
xlim([-30 30])
title('AirfoilFFAW3Cm');
