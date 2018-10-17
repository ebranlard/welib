%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied


%% Main Parameters
global Algo
dt=0.1; % time step [s]
tmax=50; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
%Controller.fpitch=@(x) (x > 10 & x<20)*2;
%Generator.fMoment=@(x) 0;
V0=10;
fExtremeOperatingGust=@(t) [0;0;V0]+[0;0; (30<t&t<40).*(-0.37*V0*sin(3*pi*(t-30)/10)).*(1-cos(2*pi*(t-30)/10))    ];
fCoherentGust=@(t)  [0;0;V0]+[0;0;(30<t&t<40).*(0.5*V0*( 1-cos(pi*(t-30)/10) ) )+(40<=t)*10 ];
%plot(0:0.1:40,fCoherentGust(0:0.1:40))





%% no pitch control
InitSystemProperties;
flag='EOG';Aero.Wind.fV0=fExtremeOperatingGust;
% flag='ECG';Aero.Wind.fV0=fCoherentGust;
% Initial condition (equilibium at 10 m/s)
x0=[ 0.077 ; 0 ;  0 ]; 
v0=[ 0 ; 2.3; 0 ];    
a0=[ 0 ; 0.0 ;  0 ];  
Runge
Thrustnp=Thrust;
Powernp=Power;
MGnp=MG;
Torquenp=Torque;
omeganp=v(2,:);

%% pitch
InitSystemProperties;
flag='EOG';Aero.Wind.fV0=fExtremeOperatingGust;
% flag='ECG';Aero.Wind.fV0=fCoherentGust;
% Initial condition (equilibium at 10 m/s)
x0=[ 0.077 ; 0 ;  0 ]; 
v0=[ 0 ; 2.3; 0 ];    
a0=[ 0 ; 0.0 ;  0 ];  
Controller.fpitch=@PitchControl;
Runge
Thrustp=Thrust;
Powerp=Power;
MGp=MG;
Torquep=Torque;
omegap=v(2,:);

%% PostProcessing
XLIM=[25 45]

setFigureWidth('17')
setFigureHeight('15')

figure
subplot(3,1,1);
hold all
plot(t,Thrustnp/1000,'b--');grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
plot(t,Thrustp/1000,'k-');grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
box on
legend('Zero pitch ','Pitch Controller',2)

subplot(3,1,2);
hold all
plot(t,Powernp/1000,'b--');grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
plot(t,Powerp/1000,'k-');grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
ylim([0 6500])
box on
%legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)

subplot(3,1,3);
hold all
plot(t,HubV);grid on;xlabel('t [s]');ylabel('Hub wind speed [m/s]');xlim(XLIM);
box on
%legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)
title([flag 'GustResponse'])
