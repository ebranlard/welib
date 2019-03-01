function [a aprime phi c]=getOptimizedParametersClCd(Cl,Cd,lambda,B)
global Rotor;
R=Rotor.R;
r=Rotor.r(:);
rhub=min(r);
lambda_r=lambda*r/R;
Cl_opt=zeros(Rotor.ne,1);
Cd_opt=zeros(Rotor.ne,1);
a=zeros(Rotor.ne,1);
for e = 1:Rotor.ne
    Cl_opt(e,1) = Cl;
    Cd_opt(e,1) = Cd;
    a(e,1)=fzero(@(aa) getRelationAX(aa,lambda_r(e)),0.3);
end
aprime=(1-3*a)./(4*a-1);
phi=atan( (1-a)./((1+aprime).*lambda_r(:)) )*180/pi;  %deg
Cn=Cl_opt.*cosd(phi)+Cd_opt.*sind(phi);

f=B/2*(R-r)./(r.*sind(phi));
Fhub=2/pi*acos(exp(-B/2*(r-rhub)./(rhub*sind(phi))));
Ftip=2/pi*acos(exp(-f));
F=Ftip.*Fhub;


c=(8*pi*R.*F.*a.*lambda_r.*(sind(phi)).^2)./((1-a)*B*lambda.*Cn);
end