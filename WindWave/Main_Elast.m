%%  Aeroelast Modelling
disp('Aeroelastic modelling')
clear all
include_globals;
% -------------- Part1
MFloater   = 7466330    ;  % kg
zCMFloater = -89.92     ;  % m
IFloaterCM = 4229230000 ;  % kg m2
MTower     = 249718     ;  % kg
zCMTower   = 43.4       ;  % m
ITowerCM   = 1.217*10^8 ;  % kg m2
MNacelle   = 240 * 10^3 ;  % kg
MRotor     = 110 * 10^3 ;  % kg
zHub       = 90         ;  % m
DSpar = 9.4    ;  % m
zBot  = -113.4 ;  % m zBottom
Cm    = 1.0    ; 
CD    = 0.60   ; 
rhow  = 1025   ;  % kg/m3     % rho water!!!
g     = 9.81   ;  % m/s2
kMoor = 41180  ;  % N/m
zMoor = -70    ;  % m
zBuoy=zBot/2;
%------------- Parameters Monopile
h=200; % m, sea bed depth in meter
%------------- Algo parameters
bWheeler=0;
bFloating=1;
bVzTo0=1;
nz=50;
bNoForces=0;
bNoWind=0; 

% ------------- mass stiffness and damping−matrix
M0 =[ 8.0660e+006   -629.0346e+006
-629.0346e+006    68.0261e+009];
A =[ 8.0665e+006  -457.3678e+006
-457.3678e+006    34.5770e+009];
M=M0+A;
C =[ 41.1800e+003    -2.8826e+006
-2.8826e+006     1.8899e+009];
damping=1.*10^5;
B= [ damping 0 ; 0 0 ];

% Time
% time
dt=0.1;
vt=0:dt:900;
istart=whichvalue(vt,300);
Iselect=istart:length(vt);
load('StochData.mat'); % contains:'vUw_Kaimal','vf_Jonswap','a_Jonswap','vphases_Jonswap','vk'

% 
vUw=10+vt*0; vf=1/10; vphases=0; vA=0;
vk=fgetDispersion(vf,h,g);
% gamma_CT=1; qInit=zeros(5,1); 
qInit=zeros(4,1); % initial condition 
[t10,q10]=ode45('fqPrime',vt,qInit);
vUw=16+vt*0; vf=1/10; vphases=0; vA=0;
vk=fgetDispersion(vf,h,g);
[t16,q16]=ode45('fqPrime',vt,qInit);
fPlotCompare(t10,q10,t16,q16,'10 m/s','16m/s','TwoWS')

%% stochastic wind stochastic wave
vUw=vUw_Kaimal+10;
vf=vf_Jonswap;
vphases=vphases_Jonswap;
vA=vA_Jonswap;
vk=fgetDispersion(vf,h,g);
% Initial condition: rest
qInit=zeros(4,1);
[tswsw,qswsw]=ode45('fqPrime',vt,qInit);

fPlotCompare(t10,q10,tswsw,qswsw,'10 m/s','10m/s SS','SS')
% [ S,f ] = fSpectrum( tswsw(Iselect),qswsw(Iselect,1),1,{' stoc. wave',' stoc. wind','x'});
% [ S,f ] = fSpectrum( tswsw(Iselect),qswsw(Iselect,2),1,{' stoc. wave',' stoc. wind','\theta'} );
