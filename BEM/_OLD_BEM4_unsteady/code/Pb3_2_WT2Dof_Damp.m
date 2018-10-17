%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='2DOF';

%% Main Parameters
global Algo
dt=0.1; % time step [s]
tmax=10; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=0;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=1;
%  Initial condition (equilibium at 10 m/s)
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

Runge
Thrustnp=Thrust;
Powernp=Power;
MGnp=MG;
Torquenp=Torque;
omeganp=v(2,:);


%% PostProcessing
XLIM=[0 max(t)]
PostPro
% 
% setFigureWidth('17')
% setFigureHeight('15')
% 
% figure
% subplot(3,1,1);
% hold all
% plot(t,Thrustnp/1000,'b--');grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
% plot(t,Thrustp/1000,'k-');grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
% box on
% legend('Zero pitch ','Pitch Controller',2)
% 
% subplot(3,1,2);
% hold all
% plot(t,Powernp/1000,'b--');grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
% plot(t,Powerp/1000,'k-');grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
% ylim([0 6500])
% box on
% %legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)
% 
% subplot(3,1,3);
% hold all
% plot(t,HubV);grid on;xlabel('t [s]');ylabel('Hub wind speed [m/s]');xlim(XLIM);
% box on
% %legend('No Generator','Generator','Pitch step','Pitch Controller','Orientation','horizontal',2)
% title([flag 'GustResponse'])
