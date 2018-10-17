%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='';
%% Main Parameters
global Algo
dt=0.01; % time step [s]
tmax=50; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
Algo.AutoStop=1;


%Vs=3.5:0.5:25;
Vs=[2:6 8:2:24];
Cases=5;

omegasAll=zeros(Cases,length(Vs));
PowersAll=zeros(Cases,length(Vs));
PowerNorm=zeros(Cases,length(Vs));
ForceNorm=zeros(Cases,length(Vs));
ThrustsAll=zeros(Cases,length(Vs));
lambdasAll=zeros(Cases,length(Vs));

A=(pi*Rotor.R^2);
Last=@(x) x(end);
LastNonZero=@(x) Last(x( x>0 | x<0  ))  ;
%% Constant omega static aerodata
jj=1;
Algo.ConstantOmega=1;
Algo.dynastall=0;

for ii=1:length(Vs)
    disp(['Current Wind speed :' num2str(Vs(ii))])
    Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
    Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
    Aero.Wind.Model='PowerLaw'; %;
    
    % Initial condition (equilibium at 10 m/s)
    x0=[ 0.0077*Vs(ii) ; 0 ;  0 ];
    v0=[ 0 ; 2.295; 0 ];
    a0=[ 0 ; 0.0 ;  0 ];
    
    % System resolution with Runge Kutta nystrom scheme
    Runge
    % Storing final equilibrium values
    PowersAll(jj,ii)=LastNonZero(Power);
    ThrustsAll(jj,ii)=LastNonZero(Thrust);
    omegasAll(jj,ii)=LastNonZero(v(2,:));
end
A=(pi*Rotor.R^2);
lambdasAll(jj,:)=omegasAll(jj,:)*Rotor.R./Vs;
PowerNorm(jj,:)=(0.5*Environment.rho*A*Vs.^3);
ForceNorm(jj,:)=(0.5*Environment.rho*A*Vs.^2);
%% Constant omega dynamic aerodata
jj=2;
InitSystemProperties;
Algo.ConstantOmega=1;
Algo.dynastall=1;

for ii=1:length(Vs)
    disp(['Current Wind speed :' num2str(Vs(ii))])
    Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
    Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
    Aero.Wind.Model='PowerLaw'; %;
    
    % Initial condition (equilibium at 10 m/s)
    x0=[ 0.0077*Vs(ii) ; 0 ;  0 ];
    v0=[ 0 ; 0.04*Vs(ii)+1.94; 0 ];
    v0=[ 0 ; 2.295; 0 ];
    a0=[ 0 ; 0.0 ;  0 ];
    
    % System resolution with Runge Kutta nystrom scheme
    Runge
    % Storing final equilibrium values
    PowersAll(jj,ii)=LastNonZero(Power);
    ThrustsAll(jj,ii)=LastNonZero(Thrust);
    omegasAll(jj,ii)=LastNonZero(v(2,:));
end
A=(pi*Rotor.R^2);
lambdasAll(jj,:)=omegasAll(jj,:)*Rotor.R./Vs;
PowerNorm(jj,:)=(0.5*Environment.rho*A*Vs.^3);
ForceNorm(jj,:)=(0.5*Environment.rho*A*Vs.^2);
%% Different omega static aerodata
jj=3;
InitSystemProperties;
Algo=rmfield(Algo,'ConstantOmega');
Algo.dynastall=0;

for ii=1:length(Vs)
    disp(['Current Wind speed :' num2str(Vs(ii))])
    Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
    Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
    Aero.Wind.Model='PowerLaw'; %;
    
    % Initial condition (equilibium at 10 m/s)
    x0=[ 0.0077*Vs(ii) ; 0 ;  0 ];
    v0=[ 0 ; 0.04*Vs(ii)+1.94; 0 ];
    v0=[ 0 ; 2.295; 0 ];
    a0=[ 0 ; 0.0 ;  0 ];
    
    % System resolution with Runge Kutta nystrom scheme
    Runge
    % Storing final equilibrium values
    PowersAll(jj,ii)=LastNonZero(Power);
    ThrustsAll(jj,ii)=LastNonZero(Thrust);
    omegasAll(jj,ii)=LastNonZero(v(2,:));
end
A=(pi*Rotor.R^2);
lambdasAll(jj,:)=omegasAll(jj,:)*Rotor.R./Vs;
PowerNorm(jj,:)=(0.5*Environment.rho*A*Vs.^3);
ForceNorm(jj,:)=(0.5*Environment.rho*A*Vs.^2);
%% Different omega dynamic aerodata
jj=4;
InitSystemProperties;
%Algo=rmfield(Algo,'ConstantOmega');
Algo.dynastall=1;

for ii=1:length(Vs)
    disp(['Current Wind speed :' num2str(Vs(ii))])
    Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
    Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
    Aero.Wind.Model='PowerLaw'; %;
    
    % Initial condition (equilibium at 10 m/s)
    x0=[ 0.0077*Vs(ii) ; 0 ;  0 ];
    v0=[ 0 ; 0.04*Vs(ii)+1.94; 0 ];
    v0=[ 0 ; 2.295; 0 ];
    a0=[ 0 ; 0.0 ;  0 ];
    
    % System resolution with Runge Kutta nystrom scheme
    Runge
    % Storing final equilibrium values
    PowersAll(jj,ii)=LastNonZero(Power);
    ThrustsAll(jj,ii)=LastNonZero(Thrust);
    omegasAll(jj,ii)=LastNonZero(v(2,:));
end
A=(pi*Rotor.R^2);
lambdasAll(jj,:)=omegasAll(jj,:)*Rotor.R./Vs;
PowerNorm(jj,:)=(0.5*Environment.rho*A*Vs.^3);
ForceNorm(jj,:)=(0.5*Environment.rho*A*Vs.^2);

%% Pitch Controller
jj=5;
InitSystemProperties;
%Algo=rmfield(Algo,'ConstantOmega');
Algo.dynastall=1;
%controller
Controller.fpitch=@(t) PitchControl(t);

for ii=1:length(Vs)
    disp(['Current Wind speed :' num2str(Vs(ii))])
    Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
    Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
    Aero.Wind.Model='PowerLaw'; %;
    
    % Initial condition (equilibium at 10 m/s)
    x0=[ 0.0077*Vs(ii) ; 0 ;  0 ];
    v0=[ 0 ; 0.04*Vs(ii)+1.94; 0 ];
    v0=[ 0 ; 2.295; 0 ];
    a0=[ 0 ; 0.0 ;  0 ];
    
    % System resolution with Runge Kutta nystrom scheme
    Runge
    % Storing final equilibrium values
    PowersAll(jj,ii)=LastNonZero(Power);
    ThrustsAll(jj,ii)=LastNonZero(Thrust);
    omegasAll(jj,ii)=LastNonZero(v(2,:));
end
A=(pi*31^2);
lambdasAll(jj,:)=omegasAll(jj,:)*Rotor.R./Vs;
PowerNorm(jj,:)=(0.5*Environment.rho*A*Vs.^3);
ForceNorm(jj,:)=(0.5*Environment.rho*A*Vs.^2);

%% Post Pro
%load(PowerFinal.mat)
save('PowerFinal2')

load('../data/PowerCurve')

%load('PowerFinal')


%%
setFigureWidth('24')
setFigureHeight('12')
figure
hold all
%set(gcf,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1], 'DefaultAxesLineStyleOrder','-|--|:')
col1=0.55*[1 1 1];
col2=0.2*[1 1 1];
plot(Vs,PowersAll(1,:)./1000,'g--+','Color',col1);
plot(Vs,PowersAll(2,:)./1000,'g--o','Color',col1);
plot(Vs,PowersAll(3,:)./1000,'r-+','Color',col2);
plot(Vs,PowersAll(4,:)./1000,'r-o','Color',col2);
plot(Vs,PowersAll(5,:)./1000,'k-o','LineWidth',2);
plot(PowerCurve(:,1),PowerCurve(:,2)/1000,'b-','LineWidth',2)
grid on;xlabel('Wind speed [m/s]');ylabel('Power [kW]');
title('PowerPowerWS');
legend('\omega cst - Static C_l','\omega cst - Dynamic C_l','\omega var - Static C_l','\omega var - Dynamic C_l','\omega var - Pitch control','Guaranteed P_E','Location','EastOutside')
ylim([0 7000])
xlim([0 24])
box on

%%
%         dlambda=dr*lambda/R;
%         CPlambda_theory(k,j)=8/lambda^2*sum(aprime.*(1-a).*x.^3.*dlambda);

load('../data/BladeDesign/IdealNoDrag')
load('../data/BladeDesign/IdealDragInfB')
load('../data/BladeDesign/IdealDragB3')

%%
setFigureWidth('0')
setFigureHeight('0')
A=(pi*31^2);
PowerNorm=(0.5*Environment.rho*A*Vs.^3);
PN2=(0.5*Environment.rho*A*PowerCurve(:,1).^3);
figure
hold all

plot(Vs,PowersAll(1,:)./PowerNorm(1,:),'g--+','Color',col1);
plot(Vs,PowersAll(2,:)./PowerNorm(1,:),'g--o','Color',col1);
plot(Vs,PowersAll(3,:)./PowerNorm(1,:),'r-+','Color',col2);
plot(Vs,PowersAll(4,:)./PowerNorm(1,:),'r-o','Color',col2);
plot(Vs,PowersAll(5,:)./PowerNorm(1,:),'k-o','LineWidth',2);
plot(PowerCurve(:,1),PowerCurve(:,2)./PN2,'b-','LineWidth',2)
grid on;xlabel('Wind speed [m/s]');ylabel('Power coefficient C_p [.]');
ylim([0 0.6])
title('PowerCpWS');
%legend('\omega cst - Static C_l','\omega cst - Dynamic C_l','\omega var - Static C_l','\omega var - Dynamic C_l','\omega var - Pitch control',0)
%%
lambdasAll

figure
hold on
plot(IdealNoDrag(:,1),IdealNoDrag(:,2),'b')
plot(IdealDragInfB(:,1),IdealDragInfB(:,2),'b--')
plot(IdealDragB3(:,1),IdealDragB3(:,2),'b-.')

plot(lambdasAll(1,:),PowersAll(1,:)./PowerNorm(1,:),'g--+','Color',col1);
plot(lambdasAll(2,:),PowersAll(2,:)./PowerNorm(1,:),'g--o','Color',col1);
plot(lambdasAll(3,:),PowersAll(3,:)./PowerNorm(1,:),'r-+','Color',col2);
plot(lambdasAll(4,:),PowersAll(4,:)./PowerNorm(1,:),'r-o','Color',col2);
plot(lambdasAll(5,:),PowersAll(5,:)./PowerNorm(1,:),'k-o','LineWidth',2);
ylim([0 0.6])
xlim([0 20])
grid on;xlabel('Tip speed ratio \lambda_R [.]');ylabel('Power coefficient C_p [.]');
title('PowerCpLambdas');
legend('Without drag','With drag and B=\infty','With drag and B=3','Location','South')
%%
figure
hold on
plot(Vs,omegasAll(1,:)*60/(2*pi),'g--+','Color',col1);
plot(Vs,omegasAll(2,:)*60/(2*pi),'g--o','Color',col1);
plot(Vs,omegasAll(3,:)*60/(2*pi),'r-+','Color',col2);
plot(Vs,omegasAll(4,:)*60/(2*pi),'r-o','Color',col2);
plot(Vs,omegasAll(5,:)*60/(2*pi),'k-o','LineWidth',2);
grid on;xlabel('Wind speed [m/s]');ylabel('Rotor speed [RPM]');
title('PowerOmegaWS');
legend('\omega cst - Static C_l','\omega cst - Dynamic C_l','\omega var - Static C_l','\omega var - Dynamic C_l','\omega var - Pitch control',2)
%%
figure
hold on
plot(lambdasAll(1,:),omegasAll(1,:)*60/(2*pi),'g--+','Color',col1);
plot(lambdasAll(2,:),omegasAll(2,:)*60/(2*pi),'g--o','Color',col1);
plot(lambdasAll(3,:),omegasAll(3,:)*60/(2*pi),'r-+','Color',col2);
plot(lambdasAll(4,:),omegasAll(4,:)*60/(2*pi),'r-o','Color',col2);
plot(lambdasAll(5,:),omegasAll(5,:)*60/(2*pi),'k-o','LineWidth',2);
grid on;xlabel('Lambda [.]');ylabel('Rotor speed [RPM]');
title('PowerOmegaLambda');
xlim([0 20])
legend('\omega cst - Static C_l','\omega cst - Dynamic C_l','\omega var - Static C_l','\omega var - Dynamic C_l','\omega var - Pitch control',1)
%%
figure
hold on
plot(Vs,lambdasAll(1,:),'g--+','Color',col1);
plot(Vs,lambdasAll(2,:),'g--o','Color',col1);
plot(Vs,lambdasAll(3,:),'r-+','Color',col2);
plot(Vs,lambdasAll(4,:),'r-o','Color',col2);
plot(Vs,lambdasAll(5,:),'k-o','LineWidth',2);
grid on;xlabel('Wind speed [m/s]');ylabel('Tip speed ratio \lambda_R [.]');
title('PowerLambdaWS');
legend('\omega cst - Static C_l','\omega cst - Dynamic C_l','\omega var - Static C_l','\omega var - Dynamic C_l','\omega var - Pitch control',1)


