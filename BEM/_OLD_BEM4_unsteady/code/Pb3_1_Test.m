%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving equilibitum induced velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
InitBEM

%% Parameters for the algorithm
V0=[8];
omega=2.295;           %rad/s
% omega=2.314;           %rad/s
dt=pi/180/omega*10;    %time step [s]
dt=0.05;
Tstop=30;              %simulation time  [s]
Vtime=0:dt:(Tstop);
nt=length(Vtime);
w_guess=-2.5;   %initial guess for induced velocity [m/s]
Model='PowerLaw'; figoff=600;
BigStorage=0;
% Steady parameters
tilt=0;
cone=0;
% Unsteady parameters
Vyaw_of_t=0;
Vpitch_of_t=Vtime*0;
VV0_of_t(1,:)=[0 0 V0];
%% With yaw model
YawModel=0;
tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UnsteadyBEM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc()


%% Plotting
% Vtime=Vtime(floor(20/dt):end);
% Power=Power(floor(20/dt):end);
% Thrust=Thrust(floor(20/dt):end);

setFigurePath('../figs/')
setFigureTitle(0)

figure
plot(Vtime,Power/1000)
xlabel('Time [s]')
ylabel('Power [kW]')
title('1power')
grid on
[29 30]
figure
plot(Vtime,Thrust)
xlabel('Time [s]')
ylabel('Thrust [kN]')
title('1thrust')
grid on
[29 30]

