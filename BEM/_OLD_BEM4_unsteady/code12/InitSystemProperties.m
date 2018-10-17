%%
global Shaft
global Nacelle
global Tower
global Rotor
global Aero
global Environment
global Generator
global Controller

%% Shaft
Shaft.Lshaft=7;       % Shaft length
Shaft.rs_in2=[0; 0; -Shaft.Lshaft];
Shaft.k=4*10^7; % [Nm/rad]

%% Generator
Generator.I=7.95*10^5;  % [kg/m^2]
Generator.fMoment=@(x) 2*10^6 *(max(x,2.1)-2.1);

%% Nacelle
Nacelle.M=53000;  % [kg]
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

%% Rotor
Rotor.R=30.56;        % Blade length [m]
Rotor.rb_center_in4=[0; 0; 0];
Rotor.M=3*9000; % mass [kg]
Rotor.I=3.4*10^6;% [kg/m/m]
%Rotor.I=3*trapz([0;Rotor.r],[0;Rotor.r.^2.*Rotor.Blade.Mass])
Rotor.cone=0;  % cone angle [deg]

%Loading reference file for geometry
fid = fopen('../data/tjaereblade.dat', 'r');
    Buffer = textscan(fid, '%f %f %f %s %s');
fclose(fid);

Rotor.r=Buffer{1};     %radial positions on the blades [m]
Rotor.chord=Buffer{2}; %chord[m]
Rotor.beta=Buffer{3};  %beta in [deg]
Rotor.ne=length(Rotor.r);    %number of element [.]
Rotor.nB=3;            %number of blades

nB=Rotor.nB;
ne=Rotor.ne;

% Loading profiles
for e=1:Rotor.ne
    M(:,:,e)=load(['../data/' eval(Buffer{4}{e})]);
    Profiles.alpha(e,:)=M(:,1,e);
    Profiles.Cl(e,:)=M(:,2,e);
    Profiles.Cd(e,:)=M(:,3,e);
    Profiles.Cl_inv(e,:)=M(:,5,e);
    Profiles.Cl_fs(e,:)=M(:,6,e);
    Profiles.f_st(e,:)=M(:,4,e);
end
Rotor.Profiles=Profiles;
% Element used for the estimate of khi (at r/R=70%)
[m Rotor.e_ref_for_khi]=min(abs(Rotor.r/Rotor.R-0.7));

% loading modes and elastic properties (they have different radial position)

% Data=load('../data/TjaerborgUsedInterpolatedData.mat','-ASCII');
% Rotor.Blade.Mass=[0;Data(:,9)];
% Rotor.Blade.r=[0;Data(:,1);]
% Rotor.Blade.eigen1e=load('../data/EigenMore/eigen1e.dat','-ASCII');
% Rotor.Blade.eigen1f=load('../data/EigenMore/eigen1f.dat','-ASCII');
% Rotor.Blade.eigen2f=load('../data/EigenMore/eigen2f.dat','-ASCII');
Data=dlmread('../data/TjaereborgAirfoilData.csv',',',1,0);
Rotor.Blade.Mass=Data(:,9);
Rotor.Blade.r=Rotor.r;
Rotor.Blade.eigen1e=zeros(ne,3);
Rotor.Blade.eigen1f=zeros(ne,3);
Rotor.Blade.eigen2f=zeros(ne,3);
Rotor.Blade.eigen1e(:,2:3)=load('../data/eigen1e.dat','-ASCII');
Rotor.Blade.eigen1f(:,2:3)=load('../data/eigen1f.dat','-ASCII');
Rotor.Blade.eigen2f(:,2:3)=load('../data/eigen2f.dat','-ASCII');
Rotor.Blade.eigenF=[1.2146 2.3908 3.5659];


for e=1:ne
   [mmm Rotor.ee(e)]=min(abs(Rotor.Blade.r-Rotor.r(e)));
end

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


%% Environment
Environment.g=9.81;         % gravity [m/s^2]
Environment.rho=1.225;      %air density [kg/m^3]


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
M_tot=Rotor.M+Nacelle.M;



