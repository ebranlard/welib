%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='12';

%% Main Parameters
global Algo
dt=0.02; % time step [s]
tmax=15; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
Algo.Weight=0;
Algo.DontRotate=0;
Controller.fpitch=@PitchControl;

uz=zeros(1,length(Vs));uy=zeros(1,length(Vs));
Vs=7:30;
%%
Vs=1:6;
for ii=1:length(Vs)
    disp(['Current Wind speed :' num2str(Vs(ii))])
    InitSystemProperties
    Aero.Wind.V0=[0; 0; Vs(ii) ]; % wind speed at hub height
    Aero.Wind.fV0=@(x)[0; 0; Vs(ii)] ; % wind speed at hub height
    Aero.Wind.Model='Constant'; %;
    
    % Initial condition (equilibium at 10 m/s)
    x0=[ 0.0077*Vs(ii) ; 0 ;  0 ;zeros(9,1)];
    v0=[ 0 ; 2.295; 0 ;zeros(9,1)];
    a0=[ 0 ; 0.0 ;  0 ;zeros(9,1)];
    
    % System resolution with Runge Kutta nystrom scheme
    Runge
    b=Rotor.Blade;
    r=Rotor.r;
    idB=1;
    uz(ii)=max(x(idB*3+1,end)*b.eigen1f(:,3)+x(idB*3+2,end)*b.eigen1e(:,3)+x(idB*3+3,end)*b.eigen2f(:,3));
    uy(ii)=max(x(idB*3+1,end)*b.eigen1f(:,2)+x(idB*3+2,end)*b.eigen1e(:,2)+x(idB*3+3,end)*b.eigen2f(:,2));
end

%% 
uz(5)=0.54;
uz(6)=0.7;
uz(2)=0.25
figure(1)
plot(Vs,uz,'+-')
xlabel('Hub wind speed [m/s]')
ylabel('Tip deflection u_z [m]')
grid on
box on
title('TipdeflectWS')