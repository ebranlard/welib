%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving equilibitum induced velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
Init
InitBEM

%% Parameters for the algorithm
V0=[10];
dt=pi/180/omega*10;             %time step [s]
Tstop=2*pi/omega*20;           %simulation time  [s]
Vtime=0:dt:(Tstop);
w_guess=-2.5;   %initial guess for induced velocity [m/s]
Model='Constant'; figoff=600;
BigStorage=1;
% Steady parameters
tilt=0;
cone=0;
% Unsteady parameters
Vyaw_of_t=20;
Vpitch_of_t=0;
VV0_of_t(1,:)=[0 0 V0];
%% With yaw model
YawModel=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UnsteadyBEM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wyaw_saved=squeeze(MW(end,:,:,:));
W0yaw_saved=squeeze(MW0(end,:,:,:));

%% Without yaw model
YawModel=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UnsteadyBEM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wyaw_nomodel_saved=squeeze(MW(end,:,:,:));
W0yaw_nomodel_saved=squeeze(MW0(end,:,:,:));

%% Without yaw at all
YawModel=1;
Vyaw_of_t=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UnsteadyBEM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_saved=squeeze(MW(end,:,:,:));
W0_saved=squeeze(MW0(end,:,:,:));

save 'WEquilibria10.mat' W_saved W0_saved W0yaw_nomodel_saved Wyaw_nomodel_saved W0yaw_saved Wyaw_saved



%% Plotting
Interv=5*329:5*382;
figure(figoff)
plot(Vtime,P)

