function [WT]=fInitWTDefault()

%% Required Params
Rotor.nB=-1;     % number of blades
Rotor.R =-1;     % Rotor radius [m]
Rotor.BladeLength =-1;     % Rotor radius [m]
Rotor.rhub = 0;  % Hub  length [m]  %%% I need to decide, check WTperf
Rotor.cone= 0;   % cone angle [deg]
Rotor.ne = -1;   
%Rotor.M=-1;      % mass [kg]
%Rotor.I=-1;       %Inertia [kg/m/m]??? kg m^2
Rotor.rb_center_in4=[0; 0; 0];
Rotor.HubHeight=60;
Rotor.Omega=-1;  % Nominal Rotational speed [rad/s]



%% Environment
Environment.g=9.81;     % gravity [m/s^2]
Environment.rho=1.225;  %air density [kg/m^3]
Environment.KinVisc=15.68*10^-6;  %Knematic viscosity nu [m^2/s]
Environment.A = 6.67;   % shape factor % from 31784
Environment.k = 2.04;   % scale factor

%% Shaft
Shaft.Lshaft=0;   % Shaft length [m]
Shaft.rs_in2=[0; 0; -Shaft.Lshaft];
%Shaft.k=-1; % [Nm/rad] 4*10^7

%% Generator
Generator.I=-1;  %7.95*10^5 [kg/m^2]
Generator.fMoment=@(x) 0; %2*10^6 *(max(x,Rotor.Omega)-Rotor.Omega);

%% Nacelle
%Nacelle.M=-1;  % [kg]
Nacelle.tilt=0;   %tilt angle [deg]

%% Tower
Tower.H=60;   % height of the tower [m]
Tower.r1=2.125;
Tower.r2=2.375;
Tower.H1=Tower.H/2;
Tower.Xtower=[-Tower.r2 -Tower.r2 -Tower.r1 Tower.r1 Tower.r2 Tower.r2];
Tower.Ytower=[0 Tower.H1 Tower.H Tower.H Tower.H1 0];
Tower.rt_in1=[Tower.H; 0; 0];
%Tower.k=-1; % [N/rad] 2*10^6


%% Specifications
Spec.Cp_rated = -1; % estimation of Cp_max at mean wind speed from 31784
Spec.P_rated = -1;  % 500kW
Spec.V_rated_thumbs = gamma(1+1/Environment.k)*Environment.A+6; % rule of thumb mean WS + 6m/s
Spec.V_rated = (Spec.P_rated/(Spec.Cp_rated*0.5*Environment.rho*Rotor.R^2*pi))^(1/3);
Spec.Reynolds_inf = 0; %75000*Rotor.R;
Spec.Reynolds_sup = 0; %150000*Rotor.R;
Spec.Omega_rated=0;
Spec.TSR_rated=0;


%% Controller
Controller.fpitch=@(x) 0;
Controller.fyaw=@(x) 0;
Controller.lastpitch=0;

%% State
State.pitch=0;    %[deg]
State.yaw=0;      %[deg]
State.psi=0;      %[deg]
State.Aero.w_guess=-2;
State.Aero.last.W=0;     %ones(3,ne,nB).*w_guess;       %temporary induced velocity
State.Aero.last.W0=0;    %ones(3,ne,nB).*w_guess;       %temporary induced velocity
State.Aero.last.W_qs=0;  %ones(3,ne,nB).*w_guess; %temporary quasistatic induced velocity
State.Aero.last.W_int=0; % ones(3,ne,nB).*w_guess; %temporary intermediate induced velocity
State.Aero.last.fs=0;
State.Aero.Power=0;
State.Aero.lastPower=0;


%% Profiles
Profiles.n=-1;

%% Sources
Sources.Format='xblade';

%% Temporary calculation
% [M K Damping]=getMatrices();
% omega0=sqrt(norm(K)/norm(M));
% T=2*pi/omega0;
%M_tot=Rotor.M+Nacelle.M;
%Rotor.I=3*trapz([0;Rotor.r],[0;Rotor.r.^2.*Rotor.Blade.Mass])

WT.Name='Default';
WT.Shaft=Shaft;
WT.Nacelle=Nacelle;
WT.Tower=Tower;
WT.DefaultEnvironment=Environment;
WT.Generator=Generator;
WT.Controller=Controller;
WT.Rotor=Rotor;
WT.Profiles=Profiles;
WT.Spec=Spec;
WT.Sources=Sources;
WT.Transient=State;
