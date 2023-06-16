function [ vz_bar ExpFactor Y Y1 Y2 kappa] = fTheodorsenExpansion(l_bar,vx,K,nB,z_inf,ntheta,bWT)
% Compute expansion factor according to Theodorsen such that R/R_w = 1+ cs*ExpansionFactor where cs is the thrust coefficient in the far wake.. It can be obtained with fTheodorsenAtOperatingPoint
%--- OUTPUT
% ExpFactor: see above
% vz_bar= z/R : downstream dimension
% Y, Y1, Y2 : Intermediate variables, see Theodorsen's book
% kappa (or xi) : mass coefficient, integral of Goldstein's circulations
%--- INPUT
% l_bar: l/R  l=h/(2 pi)  h:pitch of the helix in the far wake
% vx: vector of dimensionless radius, e.g. linspace(0,1,100) on which K is known
% K : goldstein factor
% z_inf: location downstream that can be considered as infinity, typical value 100
% ntheta: Number of points used to compute the integral over theta, typical value 3000
% bWT : wind turbine (1) or propeller (0)

kappa=trapz(vx,K.*vx);
INFTY=z_inf/l_bar;

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
        y2(ix)=K(ix)/nB*sum_y1;
    end
    Y1(it)=trapz(vx,y2);
end
% Compute Y2(theta2) once and for all
Y2=cumtrapz(vtheta,Y1)-trapz(vtheta,Y1);    
Y=l_bar^3/4*trapz(vtheta,Y2);


vz_bar=vtheta*l_bar;
dY_fromInf=l_bar^3/4*cumtrapz(vtheta,-Y2); % for wind turbines we can go from 0 to infinity, no need to reverse
Y=l_bar^3/4*trapz(vtheta,Y2);


if(bWT)
    ExpFactor=(l_bar^3)/4/kappa*cumtrapz(vtheta,-Y2); % for wind turbines we can go from 0 to infinity, no need to reverse
    ExpFactor=ExpFactor;
    vz_bar=vtheta*l_bar;
else
    ExpFactor=(l_bar^3)/4/kappa*cumtrapz(vtheta(end:-1:1),Y2(end:-1:1)); % for wind turbines we can go from 0 to infinity, no need to reverse
    ExpFactor=ExpFactor;
    vz_bar=vtheta*l_bar;
end
end

