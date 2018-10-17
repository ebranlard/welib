%% Systems initialization
InitClear;
% 
% SpecFile='data/NORDTANK_Spec.dat';
% Files={'data/test/NTK500_pitch_structure.htc','data/test/data/NewBlade.pc','data/test/data/NewBlade82.ae'};
% InitSystemProperties; 
% fInitRotor('hawc',Files)

SpecFile='data/Tjaere/Tjaere_Spec.dat';
Files={'data/Tjaere/Tjaere_BladeGeometry.dat','data/Tjaere/Tjaere_BladeProfiles.dat'};
InitSystemProperties; 
fInitRotor('flex',Files)



path(path,'./f_aeroelastic')

%%



%% Main Parameters
dt=0.05; % time step [s]
tmax=5; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=10000;
Algo.TwoDOF=1;
Algo.DOF=3;
Algo.DontRotate=0;
%Controller.fpitch=1@(x) (x > 30 & x<50)*1;
Controller.fpitch=@(x)(fStep(x,30,50,0.3,0,1,6.5));
flag='Pitch';
% Initial condition
x0=[ 0.022 ; 0 ;  0 ]; 
v0=[ 0 ; 2.844; 0 ];    
a0=[ 0 ; 0.0 ;  0 ];  

%% System resolution with Runge Kutta nystrom scheme
Runge
%% PostProcessing
%XLIM=[5 25]
%PostPro
plot(t,Power/1000,'k','LineWidth',2)


%%

% Vtime=Vtime(floor(20/dt):end);
% Power=Power(floor(20/dt):end);
% Thrust=Thrust(floor(20/dt):end);
close all
setFigurePath('./')
setFigureTitle(0)

i30=whichmin(abs(t-30))-2;
i31=whichmin(abs(t-31))-5;
p=polyfit(t([i30 i31]), Power([i30 i31])/1000,1);
tfit=29:32;
fit=polyval(p,tfit);

figure
subplot(2,1,1)
plot(t,Pitch,'k','LineWidth',2)
grid on
ylim([0 1.2])
xlim([20 60])
ylabel('Blade Pitch \theta [deg]')

subplot(2,1,2)
hold on
plot(t,Power/1000,'k','LineWidth',2)
plot(t([i30 i31]),Power([i30 i31])/1000,'o','MarkerSize',4,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
plot(tfit,fit,'--','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Time [s]')
ylabel('Aerodynamic Power [kW]')
title('1power')
grid on
box on
xlim([20 60])
ylim([300 330])
title('DpDpitchFlex')

dPdpitch_0 = (min(Power([i30:(i31+100)]))-Power(i30-10))/1000
dQdpitch_0 =  dPdpitch_0/Rotor.Omega

text('Interpreter','latex','String',sprintf('$$\\frac{dP}{d\\theta} = -$$ %.2f [kW / deg]',abs(dPdpitch_0)), 'Position',[30 324], 'FontSize',14)


