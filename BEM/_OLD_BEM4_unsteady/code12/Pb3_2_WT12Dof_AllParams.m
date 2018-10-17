%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='12';

%% Main Parameters
global Algo
dt=0.04; % time step [s]
tmax=400; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
Algo.Weight=1;
Algo.DontRotate=0;

%real average wind 
% WS=load('../data/sonic1u10Hz.dat');
% WS=1*(WS-mean(WS))+10;
% windowSize=10;
% WS=filter(ones(1,windowSize)/windowSize,1,WS);
% WS=WS(windowSize:(end-windowSize-1));
% WS(WS<0.5)=0.5;
% 
% 
% 
% Aero.Wind.fV0=@(t)[0; 0; WS(mod( floor(t/dt),length(WS) )+1)  ] ; % wind
% speed at hub height
% Aero.Wind.fV0=@(t)[0; 0; 25/tmax*t+sin(t)*(25/tmax)*t*0.14*sqrt(2)] ; % wind speed at hub height
%real average wind 
 WS=load('../data/sonic1u10Hz.dat');
 WS=(WS-mean(WS));
 WS=WS/max(abs(WS));
 windowSize=10;
 WS=filter(ones(1,windowSize)/windowSize,1,WS);
 WS=WS(windowSize:(end-windowSize-1));
 Aero.Wind.fV0=@(t)[0; 0; max(3*WS(mod( floor(t/dt),length(WS) )+1)+sin(4*pi*t/tmax)*8+12,0.1)  ] ; % wind speed at hub height


% % Stochastic turbulence
BinOpen
txlim=420;
Aero.Wind.fTurb=@(X,Y,Z,t) [MU(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MV(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MW(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1) ];  

Aero.Wind.Model='ConstantPowerLawStochasticTowerEffect';

%controller
%Controller.fpitch=@(t) (t>0&t<10)*3.5*(t)+(t>=10)*35;
Controller.fpitch=@PitchControl;

% Initial condition
a0=[ -739.8815e-006;    23.0133e-003;  -104.3992e-003; zeros(9,1)]*1;
v0=[ -2.9467e-003;  2.67e+000;  30.4165e-003;zeros(9,1)]*1;
x0=[ 78.4824e-003; 36.8422e+000; -8.3656e-003; zeros(9,1)]*1;


%% System resolution with Runge Kutta nystrom scheme
Algo.Weight=0;
%Aero.Wind.Model='PowerLawConstant';
Runge
% PostPro
%% PostProcessing

XLIM=[0 tmax];
setFigureWidth('24')
setFigureHeight('20')
figure
subplot(4,1,1);plot(t,HubV);grid on;xlabel('t [s]');ylabel('Hub speed [m/s]');xlim(XLIM);ylim([0 25])
subplot(4,1,2);plot(t,Power/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);ylim([0 2300])
subplot(4,1,3);plot(t,Pitch);grid on;xlabel('t [s]');ylabel('Pitch [deg]');xlim(XLIM);ylim([0 25])
subplot(4,1,4);plot(t,v(2,:)*60/(2*pi));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed d\theta/dt [RPM]');
title(strcat('Control12','Main'));
%%
setFigureWidth('0')
setFigureHeight('0')
load('../data/PowerCurve')
figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end)*1.05,(Power(10:end)/1000)*0.99-20,'.','MarkerSize',5)
plot(PowerCurve(:,1),PowerCurve(:,2)/1000,'k')
legend('Simulated','Guaranted ',2)
ylabel('Power [kW]')
xlabel('Wind speed [m/s]')
title('Controller12PowerCurve')
grid on
xlim([0 25])
ylim([0 2300])
%%
setFigureWidth('0')
setFigureHeight('0')
load('../data/PowerCurve')
[Vb Thrustb]=bin(HubV(10:end),Thrust(10:end),50);

figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end),(Thrust(10:end)/1000),'.','MarkerSize',5)
plot(Vb,Thrustb/1000,'k-+','MarkerSize',5)
legend('Scatter','Bins ',2)
ylabel('Thrust [kN]')
xlabel('Wind speed [m/s]')
title('Controller12ThrustCurve')
grid on
xlim([0 25])
ylim([0 250])
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
title('Controller12PitchCurve')
grid on
xlim([0 25])
 
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
legend('Bins ',2)
ylabel('Rotational speed [RPM]')
xlabel('Wind speed [m/s]')
title('Controller12RPMCurve')
grid on
xlim([0 25])
%% Flap
setFigureWidth('0')
setFigureHeight('0')

[Vb flapb]=bin(HubV(10:end),mean(Flap(10:end,:)'),50);

figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end),Flap(10:end,1)'/1000,'.','MarkerSize',5)
plot(Vb,flapb/1000,'k-+','MarkerSize',5)
legend('Scatter','Bins ',2)
ylabel('Flap Moment [kN]')
xlabel('Wind speed [m/s]')
title('Controller12FlapCurve')
grid on
ylim([0 1800])
xlim([0 25])
%% Edge
setFigureWidth('0')
setFigureHeight('0')
[Vb edgeb]=bin(HubV(10:end),mean(Edge(10:end,:)'),50);

figure
hold all
%plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(HubV(10:end),mean(Edge(10:end,:)')/1000,'.','MarkerSize',5)
plot(Vb,edgeb/1000,'k-+','MarkerSize',5)
legend('Scatter','Bins ',2)
ylabel('Edgewise moment [kN]')
xlabel('Wind speed [m/s]')
title('Controller12EdgeCurve')
grid on
xlim([0 25])
ylim([0,350]) 
 %%
setFigureWidth('24')
setFigureHeight('15')
figure;hold all
subplot(3,1,1);plot(t,sin(x(2,:)));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Position sin(\theta) [rad]');
subplot(3,1,2);plot(t,v(2,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed d\theta/dt [rad/s]');
subplot(3,1,3);plot(t,a(2,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Acc.  d^2\theta/dt^2 [rad/s^2]');
title(strcat('Rotor',flag))
%legend('Displacement [rad]','Velocity [rad/s]', 'Acceleration
%[rad/s^2]','location','north')