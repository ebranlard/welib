%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied

%% Main Parameters
global Algo
dt=0.3; % time step [s]
tmax=100; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
BigStorage=0;

%Controller.fpitch=@(x) (x > 10 & x<20)*2;
%Generator.fMoment=@(x) 0;
Aero.Wind.V0=[0;0;10];
%Shaft.k=0;
%omega=1.8;% rotation velocity [rad/s]
% Initial condition
x0=[ 0 ; 0 ;  0 ]; 
v0=[ 0 ; 0; 0 ];    
a0=[ 0 ; 0 ;  0 ];  

%% System resolution with Runge Kutta nystrom scheme
Runge
%% PostProcessing
PostPro
