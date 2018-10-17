%% Parameters for the algorithm
V0=[10];
omega=Rotor.Omega;
r=Rotor.r;
R=Rotor.R;
beta=Rotor.beta;
dt=pi/180/omega*10;    %time step [s]
dt=0.05;
Tstop=40;              %simulation time  [s]
Vtime=0:dt:(Tstop);
nt=length(Vtime);
w_guess=-2.5;   %initial guess for induced velocity [m/s]
Model='Constant'; figoff=600;
BigStorage=0;
% Steady parameters
tilt=0;
cone=0;
% Unsteady parameters
Vyaw_of_t=0;
Vpitch_of_t=Vtime*0;
Vpitch_of_t(floor(30/dt):floor(48/dt)+1)= 1*(1-exp(-(0:361)/6))  ; %??in degrees
Vpitch_of_t((floor(48/dt)+2):floor(50/dt)+1)= (exp(-(0:39)/6))  ; %??in degrees
VV0_of_t(1,:)=[0 0 V0];
% VV0_of_t=zeros(length(Vtime),3);
% VV0_of_t(:,3)=linspace(0,25,length(Vtime));

setFigurePath('../45703 - Aeroelastic Design Of Wind Turbine/adwt-report2/figs/')
%setFigurePath('./')
setFigureTitle(0)



%% With yaw model
YawModel=0;
tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UnsteadyBEM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc()
% %% Without yaw model
% YawModel=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnsteadyBEM;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting
% Vtime=Vtime(floor(20/dt):end);
% Power=Power(floor(20/dt):end);
% Thrust=Thrust(floor(20/dt):end);

i30=whichmin(abs(Vtime-30))-2;
i31=whichmin(abs(Vtime-31))-5;
p=polyfit(Vtime([i30 i31]), Power([i30 i31]),1);
tfit=29:32;
fit=polyval(p,tfit);

figure
subplot(2,1,1)
plot(Vtime,Vpitch_of_t,'k','LineWidth',2)
grid on
ylim([0 1.2])
xlim([20 60])
ylabel('Blade Pitch \theta [deg]')

subplot(2,1,2)
hold on
plot(Vtime,Power,'k','LineWidth',2)
plot(Vtime([i30 i31]),Power([i30 i31]),'o','MarkerSize',4,'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
plot(tfit,fit,'--','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)
xlabel('Time [s]')
ylabel('Aerodynamic Power [kW]')
title('1power')
grid on
box on
xlim([20 60])
ylim([300 330])
title('DpDpitchUnsteadyBEM')

dPdpitch_0 = min(Power([i30:(i31+100)]))-Power(i30-10)
dQdpitch_0 =  dPdpitch_0/Rotor.Omega

text('Interpreter','latex','String',sprintf('$$\\frac{dP}{d\\theta} = -$$ %.2f [kW / deg]',abs(dPdpitch_0)), 'Position',[30 324], 'FontSize',14)


