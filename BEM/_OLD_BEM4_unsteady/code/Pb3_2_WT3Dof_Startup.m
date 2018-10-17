%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='StartupHelped';
flag='StartupNoGen';
%% Main Parameters
global Algo
dt=0.1; % time step [s]
tmax=40; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
BigStorage=0;

Aero.Wind.fV0=@(x)[0; 0; 25] ; % wind speed at hub height
%Controller.fpitch=@(x) (x > 10 & x<20)*2;
% Generator.fMoment=@(x) 0;
%Shaft.k=0;
%omega=1.8;% rotation velocity [rad/s]
% Initial condition
% x0=[ 0 ; 0 ;  0 ]; 
% v0=[ 0 ; 0.56 ; 0 ];    
% a0=[ 0 ; 0.0220 ;  0 ];  
% % Initial condition
x0=[ 0 ; 0 ;  0 ]; 
v0=[ 0 ; 0 ; 0 ];    
a0=[ 0 ; 0 ;  0 ]; 
% a0=[
%   -739.8815e-006
%     23.0133e-003
%   -104.3992e-003];
% v0=[
%     -2.9467e-003
%      2.2997e+000
%     30.4165e-003];
% x0=[
%     78.4824e-003
%     36.8422e+000
%     -8.3656e-003];

% Controller.yaw=10;      %[deg]
% Controller.fyaw=@(x) 10;

%% System resolution with Runge Kutta nystrom scheme
Runge
%% PostProcessing
XLIM=[0 tmax];
PostPro
