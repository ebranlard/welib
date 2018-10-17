%%
global Shaft
global Nacelle
global Tower
global Aero
global Environment
global Generator
global Controller
global BEMParam
global Algo


%% Environment
Environment.g=9.81;         % gravity [m/s^2]
Environment.rho=1.225;      %air density [kg/m^3]
Environment.A = 6.67; % shape factor % from 31784
Environment.k = 2.04; % scale factor

%% Shaft
Shaft.Lshaft=7;       % Shaft length
Shaft.rs_in2=[0; 0; -Shaft.Lshaft];
Shaft.k=4*10^7; % [Nm/rad]

%% Generator
Generator.I=45791;  % [kg/m^2]
Generator.fMoment=@(x) 2*10^6 *(max(x,2.844)-2.844);

%% Nacelle
Nacelle.M=15400;  % [kg]
Nacelle.tilt=0;   %tilt angle [deg]

%% Tower
Tower.H=60;   % height of the tower [m]
Tower.r1=2.125;
Tower.r2=2.375;
Tower.H1=28;
Tower.Xtower=[-Tower.r2 -Tower.r2 -Tower.r1 Tower.r1 Tower.r2 Tower.r2];
Tower.Ytower=[0 Tower.H1 Tower.H Tower.H Tower.H1 0];
Tower.rt_in1=[Tower.H; 0; 0];
Tower.k=2*10^6; % [N/rad]


%% Aero
w_guess=-2;
Aero.last.W=ones(3,nB,ne).*w_guess;       %temporary induced velocity
Aero.last.W0=ones(3,nB,ne).*w_guess;       %temporary induced velocity
Aero.last.W_qs=ones(3,nB,ne).*w_guess; %temporary quasistatic induced velocity
Aero.last.W_int=ones(3,nB,ne).*w_guess; %temporary intermediate induced velocity
Aero.last.fs=0;
Aero.Power=0;
Aero.lastPower=0;
Aero.Wind.nu=0.2; % exponent for the power law
Aero.Wind.V0=[0; 0; 10 ]; % wind speed at hub height
Aero.Wind.fV0=@(x)[0; 0; 10] ; % wind speed at hub height
Aero.Wind.Model='Constant'; %;
Aero.Cp_rated = 0.25; % estimation of Cp_max at mean wind speed from 31784
Aero.P_rated = 500000/0.93; % 500kW
Aero.V_rated_thumbs = gamma(1+1/Environment.k)*Environment.A+6; % rule of thumb mean WS + 6m/s

Aero.V_rated = (Aero.P_rated/(Aero.Cp_rated*0.5*Environment.rho*Rotor.R^2*pi))^(1/3);
Aero.Reynolds_inf = 75000*Rotor.R;
Aero.Reynolds_sup = 150000*Rotor.R;




%% Controller
Controller.pitch=0;    %[deg]
Controller.fpitch=@(x) 0;
Controller.yaw=0;      %[deg]
Controller.fyaw=@(x) 0;
Controller.lastpitch=0;



%% Temporary calculation
% [M K Damping]=getMatrices();
% omega0=sqrt(norm(K)/norm(M));
% T=2*pi/omega0;
%M_tot=Rotor.M+Nacelle.M;


%% BEM param - default
BEMParam.nbIt=200;
BEMParam.alphaCrit=0.05;
BEMParam.relaxation=0.3;
BEMParam.correction='Spera';
BEMParam.BigStorage=0;

%% Algo Param - default
Algo.Format='HAWC';


clear('w_guess')


