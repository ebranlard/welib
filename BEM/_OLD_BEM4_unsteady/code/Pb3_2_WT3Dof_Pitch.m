%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied


%% Main Parameters
global Algo
dt=0.1; % time step [s]
tmax=90; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=0;
Algo.TwoDOF=1;
Algo.DontRotate=0;
Controller.fpitch=@(x) (x > 50 & x<70)*2;
%Generator.fMoment=@(x) 0;
%plot(0:0.1:40,fCoherentGust(0:0.1:40))
flag='Pitch';
% Initial condition
x0=[ 0.077 ; 0 ;  0 ]; 
v0=[ 0 ; 2.3; 0 ];    
a0=[ 0 ; 0.0 ;  0 ];  

%% System resolution with Runge Kutta nystrom scheme
Runge
%% PostProcessing
%XLIM=[5 25]
PostPro
