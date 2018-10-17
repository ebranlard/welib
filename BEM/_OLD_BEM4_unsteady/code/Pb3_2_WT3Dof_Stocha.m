%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='Stocha';
BinOpen
%% Main Parameters
global Algo
dt=0.05; % time step [s]
tmax=45; % simulation time [s]
% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1.1*10^6;
Algo.TwoDOF=0;

%Aero.Wind.Model='PowerLawStochasticTowerEffect'; %;
%Aero.Wind.Model='PowerLaw'; %;
% Stochastic wind loading

Aero.Wind.fTurb=@(X,Y,Z,t) [MU(floor(t*N1/50)+1, floor((Y+31)*N2/70)+1, floor((Z-Tower.H+31)*N2/70)+1)
    MV(floor(t*N1/50)+1, floor((Y+31)*N2/70)+1, floor((Z-Tower.H+31)*N2/70)+1)
    MW(floor(t*N1/50)+1, floor((Y+31)*N2/70)+1, floor((Z-Tower.H+31)*N2/70)+1) ];  

%%
%fTurb=@(X,Y,Z,t) MU(floor(t*N1/100)+1, floor((Y+31)*N2/70)+1, floor((Z-Tower.H+31)*N2/70))  ;
% figure(1)
% for i=60+(-30:5:30)
%     hold all
%     plot(0:99,fTurb(45,0,i,0:99))
% end
%%

Aero.Wind.fV0=@(x)[0; 0; 10] ; % wind speed at hub height
% Initial condition
x0=[ 0.089 ; 0 ; -0.010 ]; 
v0=[ -0.030 ; 2.317; -0.000387 ];    
a0=[ -0.234 ; -0.0027 ; 0.000718 ];  

%% System resolution with Runge Kutta nystrom scheme
Aero.Wind.Model='Constant'; %;
Runge
Pc=Power;
Tc=v(2,:);

Aero.Wind.Model='PowerLaw'; %;
Runge
Ppl=Power;
Tpl=v(2,:);

Aero.Wind.Model='TowerEffect'; %;
Runge
Pt=Power;
Tt=v(2,:);

Aero.Wind.Model='StochasticConstant'; %;
Runge
Ps=Power;
Ts=v(2,:);


Aero.Wind.Model='StochasticPowerLaw'; %;
Runge
Pspl=Power;
Tspl=v(2,:);


Aero.Wind.Model='StochasticPowerLawTowerEffect'; %;
Runge
Psplt=Power;
Tsplt=v(2,:);


%save('StochasticSmall.mat','t','Pc','Ppl','Pt','Ps','Pspl','Psplt','Tc','Tpl','Tt','Ts','Tspl','Tsplt')

%% PostProcessing
load('StochasticSmall')
XLIM=[25 45];
flag='Stocha';

%%
setFigureWidth('0')
setFigureHeight('0')
figure
hold all
plot(t,Pc/1000);
plot(t,Ppl/1000);
plot(t,Pt/1000);
legend('Constant','Shear','Tower Effect',1)
grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
ylim([700 1100])
title(strcat('Power2',flag));
%%
figure
%%hold all
set(gcf,'DefaultAxesColorOrder',meshgrid(linspace(0, 0.45,2),ones(3,1))',...
    'DefaultAxesLineStyleOrder','-')

hold all
plot(t,Ps/1000,'b');
plot(t,Pspl/1000);
plot(t,Psplt/1000,'r');
ylim([700 1100])
legend('Stochastic','Stochastic + Shear','Stochastic + Shear + Tower',1)
%legend('Constant','Shear','Tower Effect','Stochastic','Stochastic + Shear','Stochastic + Shear + Tower',0)
grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
title(strcat('Power',flag));

%%
figure;
hold all
plot(t,Tc*60/(2*pi));
plot(t,Tpl*60/(2*pi));
plot(t,Tt*60/(2*pi));
plot(t,Ts*60/(2*pi));
plot(t,Tspl*60/(2*pi));
plot(t,Tsplt*60/(2*pi));

legend('Constant','Shear','Tower Effect','Stochastic','Stochastic + Shear','Stochastic + Shear + Tower',0)
xlim(XLIM);grid on;xlabel('t [s]');ylabel('Rotor speed [RPM]');
title(strcat('Theta',flag))
%legend('Displacement [rad]','Velocity [rad/s]', 'Acceleration
%[rad/s^2]','location','north')