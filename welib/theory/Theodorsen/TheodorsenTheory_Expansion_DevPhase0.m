%% README
% This script is the developping script that I wrote to understand Theodorsen's farWake Theory and implement it
% It is meant to be user friendly
%


%% Initialization
InitClear
require('OPTIMCIRC','v01_nomex')

%% Wake Expansion / Contraction 
Rw=10;
l_bar=1;
nB=2;
w=1;
U0=10;
w_bar=w/U0;


nr=100;
ntheta_inner=300;
ntheta_outer=100;
ntheta_1=50; % watchout became useless

INFTY=50;
theta_2=0;
theta_1_max=10/l_bar;

vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
[ G kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar )
figure,
plot(vx,G)



vB=1:nB;
vtau=(vB-1)*(2*pi)/nB;



vdr=zeros(1,ntheta_1);
Y2=zeros(1,ntheta_outer);
vtheta_outer=linspace(theta_2,theta_1_max,ntheta_outer);
for ito=1:ntheta_outer
    theta_outer=vtheta_outer(ito);
    vtheta_inner=linspace(theta_outer,INFTY,ntheta_inner); 
    Y1=zeros(1,ntheta_inner);
    for iti=1:length(Y1)
        theta=vtheta_inner(iti);
        y2=zeros(1,length(vx));
        for ix=1:length(vx)
            x=vx(ix);
            sum_y1=0;
            for ib=1:length(vB)
                tau=vtau(ib);
                % !!!! I added a +10^-9 to remove the singularity
                y1=( theta*cos(theta+tau) -sin(theta+tau))*(1-2*x^2+l_bar^2*theta^2+ x*cos(theta+tau))/(1+x^2+l_bar^2*theta^2-2*x*cos(theta+tau)+10^-9)^(5/2);
                sum_y1=sum_y1+y1;
            end
            y2(ix)=G(ix)/nB*sum_y1;
        end
        Y1(iti)=trapz(vx,y2);
    end
    Y2(ito)=trapz(vtheta_inner,Y1);
    figure(112),  hold all,   plot(vtheta_inner,Y1)
end
%     vdr(it1)=-l_bar^3/4*1/kappa*trapz(vtheta_outer,Y2);
figure,plot(vtheta_outer,Y2);
vdr=cs*l_bar^3/4*1/kappa*cumtrapz(vtheta_outer,Y2);
Y=l_bar^3/4*trapz(vtheta_outer,Y2);
figure(113),  hold all,   plot(vtheta_outer,Y2)
% end
vtheta_1=vtheta_outer;
%
vz_bar=vtheta_1*l_bar;

% here's come the nastyness
R0=Rw*(1-vdr(end));
vdr=vdr*R0;

% load('vdr4.mat')
figure,hold all
plot(vz_bar,vdr/(R0*cs));

figure,hold all
plot(vz_bar,-(vdr-max(vdr))/(R0*cs));

%% New way
close all
clear all
Rw=10;
l_bar=(1);
nB=2;
w=2;
U0=10;
w_bar=w/U0;

nr=100;
ntheta=3000;

INFTY=100/l_bar;

vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
[ G kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar )
% Compute Y1(theta) once and for all, that's the computational expensive part
vtheta=linspace(0,INFTY,ntheta);
Y1=zeros(1,ntheta);
vB=1:nB;
vtau=(vB-1)*(2*pi)/nB;
for it=1:length(vtheta)
    theta=vtheta(it);
    y2=zeros(1,length(vx));
    for ix=1:length(vx)
        x=vx(ix);
        sum_y1=0;
        for ib=1:length(vB)
            tau=vtau(ib);
            % !!!! I added a +10^-9 to remove the singularity
            y1=( theta*cos(theta+tau) -sin(theta+tau))*(1-2*x^2+l_bar^2*theta^2+ x*cos(theta+tau))/(1+x^2+l_bar^2*theta^2-2*x*cos(theta+tau)+10^-32)^(5/2);
            sum_y1=sum_y1+y1;
        end
        y2(ix)=G(ix)/nB*sum_y1;
    end
    Y1(it)=trapz(vx,y2);
end
figure,plot(vtheta,Y1),title('Y1')

% Compute Y2(theta2) once and for all
Y2=cumtrapz(vtheta,Y1)-trapz(vtheta,Y1);    
figure,plot(vtheta,Y2),title('Y2')

vz_bar=vtheta*l_bar;
dY_fromInf=l_bar^3/4*cumtrapz(vtheta(end:-1:1),Y2(end:-1:1));
Y=l_bar^3/4*trapz(vtheta,Y2);
figure(113),  hold all,   plot(vtheta(end:-1:1),dY_fromInf)

figure,  hold on, grid on,   plot(vz_bar(end:-1:1),dY_fromInf./kappa)
xlim([0 10])
figure,  hold on, grid on,   plot(vz_bar(end:-1:1),dY_fromInf)
xlim([0 10])
dispatchFigs(1)
%%


%% New way adapted for wind turbines  ! to be read again since I added the WT parameters
close all
clear all
Rw=10;
l_bar=(1/3);
nB=2;
w=2;
U0=10;
w_bar=w/U0;

a=1/3;
Rw=100;
U0=10;
U=U0*(1-a);
lambda=6;
h=2*pi*Rw*(1-a)/lambda;
l=h/(2*pi);
CT=0.889;
nB=3;
Gamma=CT*pi*Rw*U0/(nB*lambda);
l_bar=h/(2*pi*Rw);
w_bar=2*a;




nr=100;
ntheta=3000;

INFTY=100/l_bar;

vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
% [ G kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar )
% 
% 

G=vx*0+Gamma*nB/(h*w_bar*U0);
kappa=trapz(vx,G.*vx);
figure,
plot(vx,G),title('K')




% Compute Y1(theta) once and for all, that's the computational expensive part
vtheta=linspace(0,INFTY,ntheta);
Y1=zeros(1,ntheta);
vB=1:nB;
vtau=(vB-1)*(2*pi)/nB;
for it=1:length(vtheta)
    theta=vtheta(it);
    y2=zeros(1,length(vx));
    for ix=1:length(vx)
        x=vx(ix);
        sum_y1=0;
        for ib=1:length(vB)
            tau=vtau(ib);
            % !!!! I added a +10^-32 to remove the singularity
            y1=( theta*cos(theta+tau) -sin(theta+tau))*(1-2*x^2+l_bar^2*theta^2+ x*cos(theta+tau))/(1+x^2+l_bar^2*theta^2-2*x*cos(theta+tau)+10^-32)^(5/2);
            sum_y1=sum_y1+y1;
        end
        y2(ix)=G(ix)/nB*sum_y1;
    end
    Y1(it)=trapz(vx,y2);
end
figure,plot(vtheta,Y1),title('Y1')

% Compute Y2(theta2) once and for all
Y2=cumtrapz(vtheta,Y1)-trapz(vtheta,Y1);    
figure,plot(vtheta,Y2),title('Y2')
%


vz_bar=vtheta*l_bar;
dY_fromInf=(l_bar^3)/4*cumtrapz(vtheta,-Y2); % for wind turbines we can go from 0 to infinity, no need to reverse
Y=l_bar^3/4*trapz(vtheta,Y2);
figure(113),  hold all,   plot(vtheta(end:-1:1),dY_fromInf)

figure,  hold on, grid on,   plot(vz_bar,dY_fromInf./kappa)
xlim([0 10])
dispatchFigs(1)



%% Using a function
close all
clear all

a=1/3;
Rw=100;
U0=10;
U=U0*(1-a);
lambda=6;
h=2*pi*Rw*(1-a)/lambda;
l=h/(2*pi);
CT=0.889;
nB=3;
Gamma=CT*pi*Rw*U0/(nB*lambda);
l_bar=h/(2*pi*Rw);
w_bar=2*a;


nr=100;
ntheta=3000;
z_inf=100;
bWT=1;

vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
% [ G kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar )

K=vx*0+Gamma*nB/(h*w_bar*U0);

[ vz_bar ExpFactor Y Y1 Y2 kappa] = fTheodorsenExpansion(l_bar,vx,K,nB,z_inf,ntheta,bWT )

figure,plot(vx,K),title('K')
figure,plot(vz_bar,Y1),title('Y1')
figure,plot(vz_bar,Y2),title('Y2')
figure,  hold on, grid on,   plot(vz_bar,ExpFactor)
xlim([0 10])
dispatchFigs(1)
%%
