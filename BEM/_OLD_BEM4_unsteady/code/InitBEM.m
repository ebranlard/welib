clc;
close all;

%% Initialization
rho=1.225;      %air density [kg/m^3]

%% Geometry 
R=30.56;        % Blade length [m]
nB=3;           % Number of Blade []
x=25;           % Radial position studied[m]
Lshaft=7;       % Shaft length
H=60;           % height of the tower [m]
omega=22/60*2*pi;% rotation velocity [rad/s]
dt=(2*pi/omega)/100;         % time step
Vpsi=0:omega*dt*180/pi:360;

VelocityParams.nu=0.2; % exponent for the power law
VelocityParams.V0=10; % wind speed at hub height


%% Translation between coordinate systems
rb_in4=[x; 0; 0];
rb_center_in4=[0; 0; 0];
rt_in1=[H; 0; 0];
rs_in2=[0; 0; -Lshaft];

%% Tower Coordinates
azim = linspace(0,2*pi,100);
Xcircl = x*cos(azim);
Ycircl = x*sin(azim)+H;
Tower.r1=2.125;
Tower.r2=2.375;
Tower.H1=28;
Tower.H=H;
Tower.Xtower=[-Tower.r2 -Tower.r2 -Tower.r1 Tower.r1 Tower.r2 Tower.r2];
Tower.Ytower=[0 Tower.H1 Tower.H Tower.H Tower.H1 0];


%Loading reference file
fid = fopen('../data/tjaereblade.dat', 'r');
Buffer = textscan(fid, '%f %f %f %s %s');
fclose(fid);

r=Buffer{1};     %radial positions on the blades [m]
chord=Buffer{2}; %chord[m]
beta=Buffer{3};  %beta in [deg]
ne=length(r);    %number of element [.]

% Loading profiles
Profiles.Cl=[];
Profiles.Cd=[];
for e=1:ne
    M(:,:,e)=load(['../data/' eval(Buffer{4}{e})]);
    Profiles.alpha(e,:)=M(:,1,e);
    Profiles.Cl(e,:)=M(:,2,e);
    Profiles.Cd(e,:)=M(:,3,e);
end
Velements=1:ne;
% Element used for the estimate of khi (at r/R=70%)
[m e_ref_for_khi]=min(abs(r/R-0.7));

% Algorithm switches
YawModel=1;
BigStorage=1;