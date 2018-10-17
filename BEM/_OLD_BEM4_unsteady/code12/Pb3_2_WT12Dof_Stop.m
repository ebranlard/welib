%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='12';

%% Main Parameters
global Algo
dt=0.02; % time step [s]
tmax=60; % simulation time [s]
Algo.tstop=20;
tpitch=10;

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^2;
Algo.TwoDOF=0;
Algo.Weight=0;
Algo.DontRotate=0;
Algo.Break=0;

tstop=Algo.tstop;


%real average wind 
WS=load('../data/sonic1u10Hz.dat');
WS=0.5*(WS-mean(WS))+20;
windowSize=10;
WS=filter(ones(1,windowSize)/windowSize,1,WS);
WS=WS(windowSize:(end-windowSize-1));
WS(WS<0.5)=0.5;
Aero.Wind.fV0=@(t)[0; 0; WS(mod( floor(t/dt),length(WS) )+1)  ] ; % wind speed at hub height
 
% Aero.Wind.fV0=@(t)[0; 0; 9];
% Stochastic turbulence
BinOpen
txlim=100;
Aero.Wind.fTurb=@(X,Y,Z,t) [MU(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MV(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MW(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1) ];  

Aero.Wind.Model='ConstantPowerLawTowerEffectStochastic';
%Aero.Wind.Model='ConstantPowerLaw';


%controller

% Controller.fpitch=@(t) (t>tstop&t<tstop+tpitch)*35/tpitch*(t-tstop)+(t>=tstop+tpitch)*35-(t>=tstop+2*tpitch)*35/tpitch*(t-tstop-2*tpitch)  +(t>=tstop+3*tpitch)*35/tpitch*(t-tstop-3*tpitch) ;

% Controller.fpitch=@(t) (t>tstop&t<tstop+tpitch)*35/tpitch*(t-tstop)+(t>=tstop+tpitch)*35;
Controller.fpitch=@(t) EmergencyPitch(t,tpitch,tstop)
%-(t>=tstop+2*tpitch)*35/tpitch*(t-tstop-2*tpitch)  +(t>=tstop+3*tpitch)*35/tpitch*(t-tstop-3*tpitch) ;

% Initial condition
a0=[ -739.8815e-006;    23.0133e-003;  -104.3992e-003; zeros(9,1)]*1;
v0=[ -2.9467e-003;  2.2997e+000;  30.4165e-003;zeros(9,1)]*1;
x0=[ 78.4824e-003; 36.8422e+000; -8.3656e-003; zeros(9,1)]*1;


%% System resolution with Runge Kutta nystrom scheme
Runge
%PostPro

%% PostProcessing
tmax=50;
t=t(1:2501);
mm=2501;
t=t(1:mm);
Power=Power(1:mm);
Pitch=Pitch(1:mm);
Torque=Torque(1:mm);
MG=MG(1:mm);
x=x(:,1:mm);
v=v(:,1:mm);
Flap=Flap(1:mm,:);
HubV=HubV(1:mm);
XLIM=[0 tmax]; 

%%
setFigureWidth('28')
setFigureHeight('35')
figure
subplot(8,1,1);plot(t,HubV);grid on;xlabel('t [s]');ylabel('Hub speed [m/s]');xlim(XLIM);ylim([0 25])
subplot(8,1,2);plot(t,Power/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);ylim([0 2300])
subplot(8,1,3);plot(t,Pitch);grid on;xlabel('t [s]');ylabel('Pitch [deg]');xlim(XLIM);ylim([0 40])
subplot(8,1,4);plot(t,Torque/1000);grid on;xlabel('t [s]');ylabel('Torque [kNm]');xlim(XLIM);
subplot(8,1,5);plot(t,v(2,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel(' d\theta/dt [rad/s]');
subplot(8,1,6);plot(t,MG/1000);grid on;xlabel('t [s]');ylabel('MG [kNm]');xlim(XLIM);
subplot(8,1,7);plot(t,mean(Flap')/1000);grid on;xlabel('t [s]');ylabel('Flap moment [kNm]');xlim(XLIM);
subplot(8,1,8);plot(t,x(4,:));grid on;xlabel('t [s]');ylabel('Flap mode 1 [arb]');xlim(XLIM);
title(strcat('Stop12','Main'));
