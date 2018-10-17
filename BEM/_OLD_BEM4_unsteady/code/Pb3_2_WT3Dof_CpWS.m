%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='';
%% Main Parameters
global Algo
dt=0.1; % time step [s]
tmax=50; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=0;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
Algo.AutoStop=1;
Algo.ConstantOmega=1;
V0=10;

%Vs=3.5:0.5:25;
Vs=4:1:25;
Vpitch=0:5:10;
%Vpitch=0;
% Vs=10
omegasStall=zeros(length(Vpitch),length(Vs));
PowersStall=zeros(length(Vpitch),length(Vs));
PowerNorm=zeros(length(Vpitch),length(Vs));
ForceNorm=zeros(length(Vpitch),length(Vs));
ThrustsStall=zeros(length(Vpitch),length(Vs));
lambdasStall=zeros(length(Vpitch),length(Vs));

A=(pi*Rotor.R^2);
Last=@(x) x(end);
LastNonZero=@(x) Last(x( x>0 | x<0  ))  ;
%%
for jj=1:length(Vpitch)
    disp(['CurrentPitch :' num2str(Vpitch(jj))])
    Controller.fpitch=@(x) Vpitch(jj);
    for ii=1:length(Vs)
        disp(['Current Wind speed :' num2str(Vs(ii))])
        Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
        Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
        Aero.Wind.Model='PowerLaw'; %;
        
        % Initial condition (equilibium at 10 m/s)
        x0=[ 0.0077*Vs(ii) ; 0 ;  0 ];
        v0=[ 0 ; 0.04*Vs(ii)+1.94; 0 ];
        v0=[ 0 ; 2.341; 0 ];
        a0=[ 0 ; 0.0 ;  0 ];
        
        % System resolution with Runge Kutta nystrom scheme
        Runge
        % Storing final equilibrium values
        PowersStall(jj,ii)=LastNonZero(Power);
        ThrustsStall(jj,ii)=LastNonZero(Thrust);
        omegasStall(jj,ii)=LastNonZero(v(2,:));
    end
    A=(pi*Rotor.R^2);
    lambdasStall(jj,:)=omegasStall(jj,:)*Rotor.R./Vs;
    PowerNorm(jj,:)=(0.5*Environment.rho*A*Vs.^3);
    ForceNorm(jj,:)=(0.5*Environment.rho*A*Vs.^2);
end
%%


% save('CpWSStall.mat','omegasStall','ThrustsStall','PowersStall','lambdasS
% tall')

%%
% load('CpWSStall.mat')
load('CPWS.mat')

%%
setFigureWidth('0')
setFigureHeight('0')

figure
plot(Vs,PowersStall./1000,'-+');
grid on;xlabel('Wind speed [m/s]');ylabel('Power [kW]');
title('PowerWS');

%%
figure
plot(Vs,PowersStall./PowerNorm,'-+');
grid on;xlabel('Wind speed [m/s]');ylabel('Power coefficient C_p [.]');
ylim([0 1])
title('CpWS');


%%
figure
plot(lambdasStall,PowersStall./PowerNorm,'-+');
ylim([0 1])
grid on;xlabel('Tip speed ratio \lambda_R [.]');ylabel('Power coefficient C_p [.]');
title('CpLambdas');

%%
figure
plot(Vs,omegasStall*60/(2*pi),'-+');
grid on;xlabel('Wind speed [m/s]');ylabel('Rotor speed [RPM]');
title('OmegaWS');
%%
figure
plot(Vs,lambdasStall,'-+');
grid on;xlabel('Wind speed [m/s]');ylabel('Tip speed ratio \lambda_R [.]');
title('LambdaWS');







%%
figure
plot(Vs,ThrustsStall./ForceNorm,'-+');
grid on;xlabel('Wind speed [m/s]');ylabel('Thrust coefficient C_T [.]');
title('CtWS');
%%
figure
plot(lambdasStall,ThrustsStall./ForceNorm,'-+');
grid on;xlabel('Tip speed ratio \lambda_R [.]');ylabel('Thrust coefficient C_T [.]');
title('CtLambdas');