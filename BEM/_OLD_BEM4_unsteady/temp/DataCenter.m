%%
M=dlmread('../data/Losses.csv',',',1)

%colormap([216 216 234; 234 232 216; 234 216 14; 216 234 13  ]/255 )
blue=[63 63 153];
bluelight=[216 216 234];
red=[153 61 113];
redlight=[234 216 226];
green=[61 153 86];
greenlight=[216 234 221];
Map=[blue; redlight; greenlight]/255
colormap(Map)
area(M(:,1),[M(:,2)';M(:,3)';M(:,4)']')
grid on
legend('Generator losses','Auxiliary power','Mechanical losses',2)
xlabel('Gross Power [kW]')
ylabel('Losses [kW]')

%%
% PE=M(:,1);
% PM=PE+M(:,2)+M(:,3);
% plot(PM,PE,PM,0.9433*PM -50)
%%
% load('../data/PowerCurve')
% PowerCurve(:,1)=PowerCurve(:,1)/1000;





%% sources of data

% martin
r3=[1.8 3.0 4.5 6 9 12 15 18 21 24 27 30];
EA3=[36 6.14 5.82 5.10 4.06 3.33 2.76 2.33 1.83 1.21 0.63 0.21];
EI13=[12000 1630 1080 623 255 129 64.8 32.4 15.2 6.04 1.82 0.32];
EI23=[12000 1725 1940 1490 905 557 349 221 131 65.7 28.1 9.5];
GIV3=[7500 362 328 207 92.8 47.7 24.7 12.9 6.23 2.57 0.84 0.18];
m3=[1700 330 389 347 283 235 196 166 172 90.3 52.6 24.2];
xe3=[0 2 54 59 63 58 51 45 41 40 47 82];
xm3=[0 2 159 165 170 158 137 121 110 102 108 136];
xs3=[0 0 11 13 18 15 15 16 17 16 14 10];
nu3=[0 5.4 0.94 1.3 1.09 0.86 0.86 0.91 0.83 0.63 0.16 -0.52];
beta3=[0 0 9 8 7 6 5 4 3 2 1 0];
AirfoilData3=[r3; EA3; EI13; EI23; GIV3 ;m3;xe3;xm3;xs3 ;nu3]'


% Stig Oye aero
r2=[3.46:1.5:30.46];
c2=[3.6 3.45 3.30 3.15 3.0 2.85 2.7 2.55 2.4 2.25 2.1 1.95 1.8 1.65 1.5 1.35 1.2 1.05 0.9];
t2=[43.54 36.25 30.58 26.53 24.1 22.55 21.13 19.85 18.7 17.69 16.81 16.07 15.46 14.92 14.38 13.84 13.30 12.76 12.22];
beta2=9:-0.5:0;
AirfoilData2=[r2; c2; beta2 ;t2]'

% stig oye dynamic
r1=[2.2 2.96 6.46 9.46 12.46 15.46 18.46 21.46 24.46 27.46 30.46];
EI11=[15000 2200 586 240 121 60.9 30.5 14.3 5.68 1.71 0.3];
EI21=[15000 2200 1400 851 524 328 208 123 61.8 26.4 8.93];
EA1=[46 7.55 4.79 3.81 3.13 2.59 2.19 1.72 1.14 0.59 0.2];
GIV1=[10000 500 207 92.8 47.7 24.7 12.9 6.23 2.57 0.84 0.18];
m1=[2100 510 390 315 250 206 173 138 94 55 25];
xe1=[0 0 59 63 58 51 45 41 40 47 82];
xm1=[0 0 165 170 158 137 121 110 102 108 136];
xs1=[0 0 13 18 15 15 16 17 16 14 10];
nu1=[0 0 1.3 1.09 0.86 0.86 0.91 0.83 0.63 0.16 -0.52];
AirfoilData1=[r1; EA1; EI11; EI21; GIV1 ;m1;xe1;xm1;xs1;nu1]'

%% a try (fail) of interpolation
%r=linspace(0.3,3,20)/3*30.56;
r=union(r1,r2);
r=union(r,r3);

beta52=interp1(r2,beta2,r);
beta53=interp1(r3,beta3,r);
c=interp1(r2,c2,r);
t=interp1(r2,t2,r);

m51=interp1(r1,m1,r);
EA51=interp1(r1,EA1,r);
EI151=interp1(r1,EI11,r);
EI251=interp1(r1,EI21,r);
GIV51=interp1(r1,GIV1,r);
xe51=interp1(r1,xe1,r);
xm51=interp1(r1,xm1,r);
xs51=interp1(r1,xs1,r);
nu51=interp1(r1,nu1,r);

m53=interp1(r3,m3,r);
EA53=interp1(r3,EA3,r);
EI153=interp1(r3,EI13,r);
EI253=interp1(r3,EI23,r);
GIV53=interp1(r3,GIV3,r);
xe53=interp1(r3,xe3,r);
xm53=interp1(r3,xm3,r);
xs53=interp1(r3,xs3,r);
nu53=interp1(r3,nu3,r);

% AirfoilData5=[r; c; t; EA51;EA53; EI151;EI153; EI251;EI253; GIV51;GIV53; m51;m53; nu51;nu53 ;xe51;xe53;xm51;xm53;xs51;xs53;beta52;beta53]'
% save('TjaerborgMegaInterpolatedData.mat','AirfoilData5','-ASCII')
%% The real values Used
r=[1.8 3 4.5 6:1.5:30];
c=interp1(r2,c2,r);
t=interp1(r2,t2,r);

m=interp1(r3,m3,r);
EA=interp1(r3,EA3,r);
EI1=interp1(r3,EI13,r);
EI2=interp1(r3,EI23,r);
GIV=interp1(r3,GIV3,r);
xe=interp1(r3,xe3,r);
xm=interp1(r3,xm3,r);
xs=interp1(r3,xs3,r);
nu=interp1(r3,nu3,r);
beta=interp1(r3,beta3,r);

% 
% r=union(r1,r2);
% r=r2
% beta=interp1(r2,beta2,r);
% c=interp1(r2,c2,r);
% t=interp1(r2,t2,r);
% 
% m=interp1(r1,m1,r);
% EA=interp1(r1,EA1,r);
% EI1=interp1(r1,EI11,r);
% EI2=interp1(r1,EI21,r);
% GIV=interp1(r1,GIV1,r);
% xe=interp1(r1,xe1,r);
% xm=interp1(r1,xm1,r);
% xs=interp1(r1,xs1,r);
% nu=interp1(r1,nu1,r);


% AirfoilData=[r; c; ;t; beta;  EA; EI1; EI2; GIV;m; xe; xm; xs; nu;beta+nu]'
% save('TjaerborgUsedInterpolatedData.mat','AirfoilData','-ASCII')
