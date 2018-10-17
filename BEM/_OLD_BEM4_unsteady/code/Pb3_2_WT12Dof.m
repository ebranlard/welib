%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='Controller';

%% Main Parameters
global Algo
dt=0.01; % time step [s]
tmax=50; % simulation time [s]

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
% BinOpen
% Aero.Wind.fTurb=@(X,Y,Z,t) [MU(mod(floor(t*N1/50),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
%     MV(mod(floor(t*N1/50),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
%     MW(mod(floor(t*N1/50),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1) ];  
% 
Aero.Wind.Model='ConstantPowerLawTowerEffect'

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
XLIM=[0 tmax];
setFigureWidth('24')
setFigureHeight('20')
figure
subplot(4,1,1);plot(t,HubV);grid on;xlabel('t [s]');ylabel('Hub speed [m/s]');xlim(XLIM);
subplot(4,1,2);plot(t,Power/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
subplot(4,1,3);plot(t,Pitch);grid on;xlabel('t [s]');ylabel('Pitch [rad]');xlim(XLIM);
subplot(4,1,4);plot(t,v(2,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed d\theta/dt [rad/s]');
title(strcat('Control',flag));
%%
load('../data/PowerCurve')
figure
hold all
plot(HubV,(Power/1000)*0.9433-50,'.','MarkerSize',5)
plot(PowerCurve(:,1),PowerCurve(:,2)/1000)
legend('Simulation')
ylabel('Power [kW]')
xlabel('Wind speed [m/s]')
title('ControllerPowerCurve')
