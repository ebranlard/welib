%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied

flag='Stoch'
%%


%% Main Parameters
global Algo
dt=0.07; % time step [s]
tmax=30; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
%Controller.fpitch=@(x) (x > 10 & x<20)*2;
%Generator.fMoment=@(x) 0;
V0=10;
%plot(0:0.1:40,fCoherentGust(0:0.1:40))

Aero.Wind.fV0=fExtremeOperatingGust;
%flag='ECG';Aero.Wind.fV0=fCoherentGust;

% Initial condition (equilibium at 10 m/s)
x0=[ 0.077 ; 0 ;  0 ]; 
v0=[ 0 ; 2.3; 0 ];    
a0=[ 0 ; 0.0 ;  0 ];  

%% System resolution with Runge Kutta nystrom scheme
Runge
%% PostProcessing
XLIM=[5 25]
PostPro
