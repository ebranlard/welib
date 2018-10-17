clear all;
clc;



%% Exercice
clear all;
clc;
r=[0.30 0.50 1.00 1.50 2.00 2.50 2.60 2.70 2.80 2.90 2.95];
r0=[0.30 0.50 1.00 1.50 2.00 2.50 2.60 2.70 2.80 2.90 2.95];
EI1=[500.00 468.75 390.62 312.50 234.38 167.19 151.56 137.50 121.88 107.81 100.00]*100;
EI2=[12000.00 11250.00 9375.00 7500.00 5625.00 4012.50 3637.50 3300.00  2925.00 2587.50 2400.00]*100;
m=[1.70 1.59 1.33 1.06 0.80 0.57 0.52 0.47 0.41 0.37 0.34];
beta=[21.32 14.66 7.69 5.17 3.31 1.38 0.90 0.42 -0.07 0.03 0.74];

% 
% r=fliplr(interp(fliplr(r),100));
% m=fliplr(interp(fliplr(m),100));
% beta=fliplr(interp(fliplr(beta),100));
% EI1=fliplr(interp(fliplr(EI1),100));
% EI2=fliplr(interp(fliplr(EI2),100));

r=linspace(0.3,3,20);
m=spline(r0,m,r);
beta=spline(r0,beta,r);
EI1=spline(r0,EI1,r);
EI2=spline(r0,EI2,r);


V0=[5 8 14 20];
py0=[2.5 15 50 70];
pz0=[90 160 260 280];


beta=beta*pi/180; % now beta is in radians
nu=0;
    
n=length(r);
R=3;

% modif
% n=200;
% r=(1:n)/n;
% EI1=ones(1,n);
% EI2=ones(1,n);
% m=ones(1,n);
% beta=zeros(1,n);
% EI1=ones(1,n)*mean(EI1);
% EI2=ones(1,n)*mean(EI2);


colrs={'k-+','b-+','r-+','g-+','m-+'};
for k=1:length(V0)
    V=V0(k);
    disp(sprintf('Calculating for wind speed %d',V));

    %Assuming constant load distribution
    py=ones(1,n)*py0(k);
    pz=ones(1,n)*pz0(k);    
    
    %getting uz and uy
    [uy uz]=fDeflectionFromLoad(py,pz,beta,nu,r,EI1,EI2);

    figure(3)
    hold on
    plot(r,uy,colrs{k})
    
    figure(4)
    hold on
    plot(r,uz,colrs{k})
end


figure(3)
xlabel('r [m]')
ylabel('uy [m]')
grid
box
legend('V = 5 m/s', 'V = 8 m/s' , 'V = 14 m/s' , 'V = 20 m/s','Location','NorthWest')

figure(4)
xlabel('r [m]')
ylabel('uz [m]')
grid
box
legend('V = 5 m/s', 'V = 8 m/s' , 'V = 14 m/s' , 'V = 20 m/s','Location','NorthWest')
