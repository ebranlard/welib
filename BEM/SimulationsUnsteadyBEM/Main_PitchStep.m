% This script used to be called Main_Pb3_1_PitchStep
% It's a mix between aeroelasticity morten and Martin
%
%% Initialization
InitClear
require('BEM','v04');
require('WTlib','v04');

sWT='SB2'; Format='xblade'; Pref=2000*1000;
% sWT='NTK500'; Format='hawc'; Pref=500*1000;
sWT='Tjaere'; Format='flex'; Pref=2200*1000;
[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
WT.Nacelle.tilt=0;   %<--------------------------------------------------
WT.Rotor.cone=0;
r0=[    6.4600 9.4600 12.4600 15.4600 18.4600 21.4600 24.4600 27.4600 28.9600 29.8600 30.5500];
% [ WT ]   = fSetRotorGrid(24,WT);
[ WT ]   = fSetRotorGrid(r0,WT);
[ Algo ] = fInitAlgo();
V0=8.7;
DeltaPitch=3.7; 
omega=WT.Spec.Omega_rated;


%% Parameters for the algorithm
Algo.bSteady=0;
[ Sim ]= fInitSim(WT, [V0 omega*60/(2*pi) 0 0]);

Algo.dt=pi/180/omega*10;    %time step [s]
Algo.dt=0.05;
Algo.tmax=90;              %simulation time  [s]  <------------------------------------should be 70 or 80
Algo.Vtime=0:Algo.dt:(Algo.tmax);

[ Wind ] = fInitWind(Sim);
Wind.Model='Constant'; figoff=600;
% Wind.V0=[0; 0; V0 ]; % wind speed at hub height
% Wind.fV0=@(x)[0; 0; V0] ; % wind speed at hub height
Algo.VV0_of_t=[0;0;V0];


% Steady parameters
% Unsteady parameters
tstart_record = 30 ; 
tstart_step   = 32 ; 
%                               t     t0       t1 tau A0 A1        k
WT.Controller.fpitch=@(x)(fStep(x,tstart_step,63.5,0.2,0,DeltaPitch,6.5));

%% BEM
Algo.BEM.Correction='Spera';
Algo.bAIDrag=0;
Algo.bUnsteadyStorage=0;
Algo.BEM.bTipLoss=1;  %< ------------------------------------------------------ should be one
Algo.BEM.bYawModel=0;
Algo.BEM.w_guess=-2.5;
Algo.BEM.bDynaWake=1;% <-------------------------------------------------- watch out
Algo.BEM.Kind='Bramwell';
% [R]= fWTSimulation('BEM',WT,Sim,Wind,Algo)
% R = fRunBEM(WT,Sim,Wind,Algo);
WT.Aero.last.W=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;      %temporary induced velocity
WT.Aero.last.W0=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;     %temporary induced velocity
WT.Aero.last.W_qs=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;   %temporary quasistatic induced velocity
WT.Aero.last.W_int=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;  %temporary intermediate induced velocity
WT.Aero.last.chi=0; % Kind of Yaw
[R WT]=fBEMunsteady(WT,Sim,Wind,Algo);


% Algo.relaxation=0;
Algo.bSteady=1;
% Algo.nbIt=1; BEM1 = fRunBEM(WT,Sim,Wind,Algo);
% Algo.nbIt=2; BEM2 = fRunBEM(WT,Sim,Wind,Algo);
% Algo.nbIt=3; BEM3 = fRunBEM(WT,Sim,Wind,Algo);
Algo.relaxation=.3;  % !!!!!!!!!!!!!!!!!!!!!!!if relaxation higher, then the torque is too big, all of this because of the tip..
Algo.bSteady=1;
Algo.nbIt=200; 
BEM = fRunBEM(WT,Sim,Wind,Algo);
BEM.Torque/1000



Power=R.Power{1}/1000;
Torque=R.Torque{1}/1000;
Thrust=R.Thrust{1}/1000;
Vtime=Algo.Vtime;




%% PLotting Torque standalone
load('TjaerPitch.mat')
% InitCell
figure
hold on
plot(TjaerPitch(:,1),TjaerPitch(:,2),'k+')
plot(Vtime-tstart_record,Torque,'k','LineWidth',2)
% plot(Vtime-tstart_record,Vtime*0+BEM.Torque/1000,'k','LineWidth',2)
% plot(Vtime([i30 i31])-tstart_record,Torque([i30 i31]),'o','MarkerSize',4,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
% plot(tfit-tstart_record,fit,'--','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Time [s]')
ylabel('Aerodynamic Torque [kNm]')
title('PitchStepTjaerBorgBEMUnsteady')
legend('Measurements','BEM unsteady')
grid on
box on
xlim([0 60])
ylim([100 400])









%% Computation of dpdpitch
i30=whichmin(abs(Algo.Vtime-30))-2;
i31=whichmin(abs(Algo.Vtime-31))-5;
p=polyfit(Algo.Vtime([i30 i31]), Power([i30 i31]),1);
tfit=29:32;
fit=polyval(p,tfit);
dPdpitch_0 = min(Torque([i30:(i31+100)]))-Torque(i30-10)
dQdpitch_0 =  dPdpitch_0/WT.Rotor.Omega;







figure
hold on
plot(Vtime,Torque,'r','LineWidth',2)
plot(Vtime,Vtime*0+BEM.Torque/1000,'k','LineWidth',2)
ylim([100 400])



%% Plotting  for the torque
InitCell
i30=whichmin(abs(Algo.Vtime-tstart_step))-2;
i31=whichmin(abs(Algo.Vtime-tstart_step-1))-5;
p=polyfit(Algo.Vtime([i30 i31]), Torque([i30 i31]),1);
tfit=29:32;
fit=polyval(p,tfit);

figure
subplot(2,1,1)
plot(Algo.Vtime,R.Pitch{1},'k','LineWidth',2)
grid on
% ylim([0 1.2])
xlim([20 Algo.tmax])
ylabel('Blade Pitch \theta [deg]')

subplot(2,1,2)
hold on
plot(Vtime,Torque,'k','LineWidth',2)
plot(Vtime([i30 i31]),Torque([i30 i31]),'o','MarkerSize',4,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
plot(tfit,fit,'--','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Time [s]')
ylabel('Aerodynamic Torque [kNm]')
title('1power')
grid on
box on
xlim([tstart_record Algo.tmax])
% ylim([300 330])
ylim([100 400])
title('DQDpitchUnsteadyBEM')
text('Interpreter','latex','String',sprintf('$$\\frac{dQ}{d\\theta} = -$$ %.2f [kW / deg]',abs(dQdpitch_0)), 'Position',[30 324], 'FontSize',14)





%% Plotting  for the power
figure
subplot(2,1,1)
plot(Algo.Vtime,R.Pitch{1},'k','LineWidth',2)
grid on
% ylim([0 1.2])
xlim([tstart_record Algo.tmax])
ylabel('Blade Pitch \theta [deg]')

subplot(2,1,2)
hold on
plot(Vtime,Power,'k','LineWidth',2)
plot(Vtime([i30 i31]),Power([i30 i31]),'o','MarkerSize',4,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
plot(tfit,fit,'--','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Time [s]')
ylabel('Aerodynamic Power [kW]')
grid on
box on
xlim([tstart_record Algo.tmax])
ylim([500 730])
title('DPDpitchUnsteadyBEM')
text('Interpreter','latex','String',sprintf('$$\\frac{dP}{d\\theta} = -$$ %.2f [kW / deg]',abs(dPdpitch_0)), 'Position',[30 654], 'FontSize',14)





%%
figure
plot(Vtime,Thrust)
xlabel('Time [s]')
ylabel('Thrust [kN]')
title('1thrust')
grid on
xlim([tstart_record Algo.tmax])
% export2pdf

%%
% 
% figure
% plot(Vtime,Torque)
% xlabel('Time [s]')
% ylabel('Torque [kNn]')
% title('1torque')
% grid on
% xlim([tstart_record Algo.tmax])
%ylim([100 400])
%%
% figure
% plot(HubV,Power,'b')
% xlabel('WS [m/s]')
% ylabel('Power [kW]')
% title('1power')
% grid on
