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

trapz([0 r1],[0 m1])
%% a try (fail) of interpolation
%r=linspace(0.3,3,20)/3*30.56;
r=[ 6.46 9.46 12.46 15.46 18.46 21.46 24.46 27.46 28.96 29.86 30.46]
     
beta=interp1(r2,beta2,r);
c=interp1(r2,c2,r);
t=interp1(r2,t2,r);

m=interp1(r1,m1,r);
EA=interp1(r1,EA1,r);
EI1=interp1(r1,EI11,r);
EI2=interp1(r1,EI21,r);
GIV=interp1(r1,GIV1,r);
xe=interp1(r1,xe1,r);
xm=interp1(r1,xm1,r);
xs=interp1(r1,xs1,r);
nu=interp1(r1,nu1,r);
trapz([0 r],[0 m])
AirfoilData=[r; c;t;beta; EA; EI1; EI2; GIV ;m;xe;xm;xs;nu;beta+nu]';
%%
%dlmwrite('../data/TjaereborgAirfoilData.csv',AirfoilData)
