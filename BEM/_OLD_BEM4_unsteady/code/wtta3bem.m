clear all;
clc;
%% our Geometry
loading;
B=3;
R=3;
rhub=0.28;
lambda=8;
alpha_d=6;  % alpha(17)=6 deg
r= [0.3 0.5 1.0 1.5 2.0 2.5 2.6 2.7 2.8 2.9 2.95 ];
%size of each element dr
l=([rhub r]+[r R])/2;l(1)=rhub;l(length(l))=R;
dr=diff(l);     
x=r/R*lambda;   % this is actually lambda_r
[a aprime phi c beta ]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x);
% phi is in rad
% beta is in deg

%% Exercice
r0= [0.3 0.5 1.0 1.5 2.0 2.5 2.6 2.7 2.8 2.9 2.95 ];
EI1=[500.00 468.75 390.62 312.50 234.38 167.19 151.56 137.50 121.88 107.81 100.00]*100;
EI2=[12000.00 11250.00 9375.00 7500.00 5625.00 4012.50 3637.50 3300.00  2925.00 2587.50 2400.00]*100;
m=[1.70 1.59 1.33 1.06 0.80 0.57 0.52 0.47 0.41 0.37 0.34];
beta=[21.32 14.66 7.69 5.17 3.31 1.38 0.90 0.42 -0.07 0.03 0.74];

% r=linspace(0.3,3,20);
% m=spline(r0,m,r);
% beta=spline(r0,beta,r);
% EI1=spline(r0,EI1,r);
% EI2=spline(r0,EI2,r);

V0=[5 8 14 20];
py0=[2.5 15 50 70];
pz0=[90 160 260 280];
omega=21.33;
nu=0;

%% Looping on wind speeds
close all
n=length(r);
colrs={'k-+','b-+','r-+','g-+','m-+'};
colrss={'k-.','b-.','r-.','g-.','m-.'};
for k=1:length(V0)
    V=V0(k);
    disp(sprintf('Calculating for wind speed %d',V));

    %% getting loads from BEM code
    lambda=omega*R/V;
    %lambda=8;
    Param.nbIt=5000;
    Param.alphaCrit=0.05;
    Param.relaxation=1;
    Param.correction='none';
    [BEM CP CT CQ] = BEMfunction(Param,Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V,rho);
    %Assuming constant load distribution
    pyy=ones(1,n)*py0(k);
    pzz=ones(1,n)*pz0(k);    
    [ trapz(r,pyy) trapz(r,pzz)]
    
    pzzz=lasts(BEM.Cn).*(0.5*rho*(BEM.Vr').^2.*c)/R;
    pyyy=lasts(BEM.Ct).*(0.5*rho*(BEM.Vr').^2.*c)/R;
    off=pi/2*0;
    py=pyyy.*cos(beta*pi/180+off)+pzzz.*sin(beta*pi/180+off);
    pz=-pyyy.*sin(beta*pi/180+off)+pzzz.*cos(beta*pi/180+off);
    
%     py=BEM.Q/B;
%     pz=BEM.T/B;
    [ trapz(r,py) trapz(r,pz)]
    
    %getting uz and uy
    [uy uz]=fDeflectionFromLoad(py,pz,beta*pi/180,nu,r,EI1,EI2);

    figure(23)
    hold on
    plot(r,uy,colrs{k})
    
    figure(24)
    hold on
    plot(r,uz,colrs{k})
    figure(25)
    hold on
    plot(r,py,colrs{k})
    plot(r,pyy,colrss{k})    
    figure(26)
    hold on
    plot(r,pz,colrs{k})
    plot(r,pzz,colrss{k})    
end


figure(23)
xlabel('r [m]')
ylabel('uy [m]')
grid
box
legend('V = 5 m/s', 'V = 8 m/s' , 'V = 14 m/s' , 'V = 20 m/s','Location','NorthWest')

figure(24)
xlabel('r [m]')
ylabel('uz [m]')
grid
box
legend('V = 5 m/s', 'V = 8 m/s' , 'V = 14 m/s' , 'V = 20 m/s','Location','NorthWest')


figure(25)
xlabel('r [m]')
ylabel('py [m]')
grid
box
legend('V = 5 m/s','', 'V = 8 m/s' , '','V = 14 m/s' ,'', 'V = 20 m/s','','Location','NorthWest')

figure(26)
xlabel('r [N/m]')
ylabel('pz [N/m]')
grid
box
legend('V = 5 m/s','', 'V = 8 m/s' , '','V = 14 m/s' ,'', 'V = 20 m/s','','Location','NorthWest')
