%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied

%% Main Parameters
flag='';
global Algo
dt=0.01; % time step [s]
tmax=100; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;

Aero.Wind.fV0=@(x)[0; 0; 10] ; % wind speed at hub height
% Initial condition
% Initial condition
x0=[ 0 ; 0 ;  0 ]; 
v0=[ 0 ; 0 ; 0 ];    
a0=[ 0 ; 0 ;  0 ]; 



%% No Generator
InitSystemProperties;
Generator.fMoment=@(x) 0;
Runge
Thrustng=Thrust;
Powerng=Power;
MGng=MG;
Torqueng=Torque;
omegang=v(2,:);

%% Generator
InitSystemProperties;
Generator.fMoment=@(x) 2*10^6 *(max(x,2.1)-2.1);
Runge
Thrustg=Thrust;
Powerg=Power;
MGg=MG;
Torqueg=Torque;
omegag=v(2,:);

%% PitchStep
InitSystemProperties;
Controller.fpitch=@(x) (x > 50 & x<70)*2;
Runge
Thrustp=Thrust;
Powerp=Power;
MGp=MG;
Torquep=Torque;
omegap=v(2,:);

%% Pitch controlled
InitSystemProperties;
Controller.fpitch=@(t) PitchControl(t);
Runge
Thrustpc=Thrust;
Powerpc=Power;
MGpc=MG;
Torquepc=Torque;
omegapc=v(2,:);

%%
load('StartupFinal')


XLIM=[0 max(t)];
setFigureWidth('26')
setFigureHeight('20')

figure
subplot(4,1,1);
hold all
plot(t,Thrustng/1000);grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
plot(t,Thrustg/1000);grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
plot(t,Thrustp/1000);grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
plot(t,Thrustpc/1000);grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
box on
legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)

subplot(4,1,2);
hold all
plot(t,Torqueng/1000);grid on;xlabel('t [s]');ylabel('Torque [kNm]');xlim(XLIM);
plot(t,Torqueg/1000);grid on;xlabel('t [s]');ylabel('Torque [kNm]');xlim(XLIM);
plot(t,Torquep/1000);grid on;xlabel('t [s]');ylabel('Torque [kNm]');xlim(XLIM);
plot(t,Torquepc/1000);grid on;xlabel('t [s]');ylabel('Torque [kNm]');xlim(XLIM);
box on
%legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)

subplot(4,1,3);
hold all
plot(t,MGng/1000);grid on;xlabel('t [s]');ylabel('MG [kNm]');xlim(XLIM);
plot(t,MGg/1000);grid on;xlabel('t [s]');ylabel('MG [kNm]');xlim(XLIM);
plot(t,MGp/1000);grid on;xlabel('t [s]');ylabel('MG [kNm]');xlim(XLIM);
plot(t,MGpc/1000);grid on;xlabel('t [s]');ylabel('MG [kNm]');xlim(XLIM);
box on
%legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)

subplot(4,1,4);
hold all
plot(t,Powerng/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
plot(t,Powerg/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
plot(t,Powerp/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
plot(t,Powerpc/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
box on
%legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)
title('StartupFinal')
