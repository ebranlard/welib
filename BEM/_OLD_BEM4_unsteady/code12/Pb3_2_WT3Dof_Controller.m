%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='Controller';

%% Main Parameters
global Algo
dt=0.01; % time step [s]
tmax=100; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
BigStorage=0;
% Algo.ConstantOmega=0;

%real average wind 
WS=load('../data/sonic1u10Hz.dat');
WS=2*(WS-mean(WS))+10;
Controller.pitch=5;
windowSize=10;
WS=filter(ones(1,windowSize)/windowSize,1,WS);
WS=WS(windowSize:(end-windowSize-1));
WS(WS<0.5)=0.5;
Aero.Wind.fV0=@(t)[0; 0; WS(mod( floor(t/dt),length(WS) )+1)  ] ; % wind speed at hub height
%  Aero.Wind.fV0=@(t)[0; 0;15 ] ; % wind speed at hub height
% Controller.fpitch=@(t) 0;

% 
% % Stochastic turbulence
BinOpen
Aero.Wind.fTurb=@(X,Y,Z,t) [MU(mod(floor(t*N1/50),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MV(mod(floor(t*N1/50),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MW(mod(floor(t*N1/50),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1) ];  

Aero.Wind.Model='ConstantPowerLawTowerEffectStochastic'

%controller
Controller.fpitch=@(t) PitchControl(t);


% Initial condition
a0=[
  -739.8815e-006
    23.0133e-003
  -104.3992e-003];
v0=[
    -2.9467e-003
     2.2997e+000
    30.4165e-003];
x0=[
    78.4824e-003
    36.8422e+000
    -8.3656e-003];

%% System resolution with Runge Kutta nystrom scheme
Runge
%% PostProcessing
%load('ControllerAll')


XLIM=[0 tmax];
setFigureWidth('24')
setFigureHeight('20')
figure
subplot(4,1,1);plot(t,HubV);grid on;xlabel('t [s]');ylabel('Hub speed [m/s]');xlim(XLIM);ylim([0 25])
subplot(4,1,2);plot(t,Power/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);ylim([0 2300])
subplot(4,1,3);plot(t,Pitch);grid on;xlabel('t [s]');ylabel('Pitch [deg]');xlim(XLIM);ylim([0 25])
subplot(4,1,4);plot(t,v(2,:)*60/(2*pi));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed d\theta/dt [RPM]');
title(strcat('Control','Main'));
%%
setFigureWidth('0')
setFigureHeight('0')
load('../data/PowerCurve')
figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end)*1.06,(Power(10:end)/1000)*0.99-50,'.','MarkerSize',5)
plot(PowerCurve(:,1),PowerCurve(:,2)/1000,'k')
legend('Simulated','Guaranted ',2)
ylabel('Power [kW]')
xlabel('Wind speed [m/s]')
title('ControllerPowerCurve')
grid on
xlim([4 22])
ylim([0 2300])
%%
setFigureWidth('0')
setFigureHeight('0')
load('../data/PowerCurve')
[Vb Thrustb]=bin(HubV(10:end),Thrust(10:end),50);

figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end)*1.06,(Thrust(10:end)/1000),'.','MarkerSize',5)
plot(Vb*1.06,Thrustb/1000,'k-+','MarkerSize',5)
legend('Scatter','Bins ',2)
ylabel('Thrust [kN]')
xlabel('Wind speed [m/s]')
title('ControllerThrustCurve')
grid on
 xlim([4 22])
%ylim[]
%%
setFigureWidth('0')
setFigureHeight('0')
[Vb Pitchb]=bin(HubV(10:end),Pitch(10:end),50);

figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end),(Pitch(10:end)),'.','MarkerSize',5)
plot(Vb,Pitchb,'k-+','MarkerSize',5)
legend('Scatter','Bins ',2)
ylabel('Pitch [deg]')
xlabel('Wind speed [m/s]')
title('ControllerPitchCurve')
grid on
 xlim([4 22])
 
%%
setFigureWidth('0')
setFigureHeight('0')
omega=v(2,:);
[Vb omegab]=bin(HubV(10:end),mod(omega(10:end),2*pi),50);

figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end),mod(omega(10:end),2*pi)*60/(2*pi),'.','MarkerSize',5)
plot(Vb,omegab*60/(2*pi),'k-+','MarkerSize',5)
legend('Scatter','Bins ',2)
ylabel('Rotational speed [RPM]')
xlabel('Wind speed [m/s]')
title('ControllerRPMCurve')
grid on
 xlim([4 22])