%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='12';

%% Main Parameters
global Algo
dt=0.05; % time step [s]
tmax=20; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
Algo.Weight=1;
Algo.DontRotate=0;


%real average wind 
 WS=load('../data/sonic1u10Hz.dat');
 WS=(WS-mean(WS));
 WS=WS/max(abs(WS));
 windowSize=10;
 WS=filter(ones(1,windowSize)/windowSize,1,WS);
 WS=WS(windowSize:(end-windowSize-1));
 Aero.Wind.fV0=@(t)[0; 0; max(3*WS(mod( floor(t/dt),length(WS) )+1)+sin(4*pi*t/tmax)*8+12,0.1)  ] ; % wind speed at hub height

% Aero.Wind.fV0=@(t)[0; 0; 25/tmax*t+sin(t)*(25/tmax)*t*0.14*sqrt(2)] ; % wind speed at hub height
%real average wind 
% Aero.Wind.fV0=@(t)[0; 0; 10  ] ; % wind speed at hub height
 
% % Stochastic turbulence
BinOpen
txlim=610;
Aero.Wind.fTurb=@(X,Y,Z,t) [MU(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MV(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1)
    MW(mod(floor(t*N1/txlim),N1)+1, mod(floor((Y+31)*N2/70),N2)+1, mod(floor((Z-Tower.H+31)*N2/70),N3)+1) ];  

Aero.Wind.Model='ConstantPowerLawStochastic'

%controller
%Controller.fpitch=@(t) (t>0&t<10)*3.5*(t)+(t>=10)*35;
Controller.fpitch=@PitchControl;

% Initial condition
a0=[ -739.8815e-006;    23.0133e-003;  -104.3992e-003; zeros(9,1)]*1;
v0=[ -2.9467e-003;  2.2997e+000;  30.4165e-003;zeros(9,1)]*1;
x0=[ 78.4824e-003; 36.8422e+000; -8.3656e-003; zeros(9,1)]*1;


%% System resolution with Runge Kutta nystrom scheme
Algo.Weight=0;
%Aero.Wind.Model='PowerLawConstant';
Runge
% PostPro

%% fd
% mm=13000;
% t=t(1:mm);
% Power=Power(1:mm);
% Pitch=Pitch(1:mm);
% x=x(:,1:mm);
% v=v(:,1:mm);
% Flap=Flap(:,1:mm);
% Edge=Edge(:,1:mm);
% HubV=HubV(1:mm);

%% PostProcessing

%%
fp=0.369;
fpp=fp*2;
fppp=fp*3;

 
%%
close all
figure(4)
setFigureTitle(1)
[f,S]=getSpectrum(t,dt,HubV);
loglog(f,sqrt(S),'Color','b');
grid on
xlabel('Frequency f [Hz]')
ylabel('Amplitude [m^2/s^2/Hz]')
title('   ')
% export2pdf



%%
close all
figure(2)
hold all
setFigureTitle(1)
[f,S]=getSpectrum(t,dt,Edge(:,2)'/1000);
line(f,sqrt(S),'Color','b');
xlim([0.01 4]);
grid on
title('   ')
xlabel('Frequency f [Hz]')
ylabel('Amplitude [kN/Hz]')
get(gca,'Position')
pos=[112.9355e-003   109e-003   878.639e-003  815.8280e-003];
ax2 = axes('Position',pos,'Color','none','XColor','k','LineWidth',0.5,'YAxisLocation','right' ,'XAxisLocation','top','XTick',[fp fpp fppp],'XTickLabel','1P|2P|3P','YTick',[]);
box on
set(ax2,'XGrid','on','GridLineStyle','-.')
xlim([0.01 4]);
% export2pdf
% title('SpectrumEdge')
%%
close all
figure(3)
hold all
setFigureTitle(1)
%[f,S]=getSpectrum(t,dt,mean(Flap(:,1)/1000));
[f,S]=getSpectrum(t,dt,Power/1000);
line(f,sqrt(S),'Color','b');
xlim([0.01 4]);
grid on
xlabel('Frequency f [Hz]')
ylabel('Amplitude [kW/Hz]')
title('   ')
get(gca,'Position')
pos=[112.9355e-003   109e-003   878.639e-003  815.8280e-003];
ax2 = axes('Position',pos,'Color','none','XColor','k','LineWidth',0.5,'YAxisLocation','right' ,'XAxisLocation','top','XTick',[fp fpp fppp],'XTickLabel','1P|2P|3P','YTick',[]);
box on
set(ax2,'XGrid','on','GridLineStyle','-.')
xlim([0.01 4]);
% export2pdf
% title('SpectrumPower')

%%
% plot(Edge(:,1))