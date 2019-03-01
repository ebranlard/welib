function [a aprime phi c beta]=getOptimizedParameters(alpha_d,lambda,B)
global Rotor Algo Profiles;
R=Rotor.R;
r=Rotor.r;
lambda_r=lambda*r/R;
Cl_opt=zeros(Rotor.ne,1);
Cd_opt=zeros(Rotor.ne,1);
a=zeros(Rotor.ne,1);
for e = 1:Rotor.ne
    if(Algo.Cl2piAlpha)
            Cl_opt=2*pi*alpha_d;
            Cd_opt=0;        
    else
        if(~isequal(Algo.Format,'flex'))
            ClCdCm= fAeroCoeff(alpha_d,Profiles,Rotor.ProfileSet(:,e),Rotor.thickness_rel(e),1);
            Cl_opt=ClCdCm(1);
            Cd_opt=ClCdCm(2);
        else
            ee=Rotor.ProfileSet(2,e);
            % Badly programmed, what if all the alphas are not the same,
            % then the use of a table is bad
            % static aerodynamic coefficients
            Cd_opt(e,1)= interp1(Profiles.alpha(:,ee) , Profiles.Cd(:,ee)  , alpha_d);
            Cl_opt(e,1)= interp1(Profiles.alpha(:,ee) , Profiles.Cl(:,ee)  , alpha_d);
        end
    end
    %alpha_data= Rotor.Profiles.alpha(Rotor.pe(e),:);
    %Cl_opt(e,1) = interp1(alpha_data,Rotor.Profiles.Cl(Rotor.pe(e),:),alpha_d);
    %Cd_opt(e,1) = interp1(alpha_data,Rotor.Profiles.Cd(Rotor.pe(e),:),alpha_d);
    
    a(e,1)=fzero(@(aa) getRelationAX(aa,lambda_r(e)),0.3);
end
aprime=(1-3*a)./(4*a-1);
phi=atan( (1-a)./((1+aprime).*lambda_r) )*180/pi;  %deg
Cn=Cl_opt.*cosd(phi)+Cd_opt.*sind(phi);

f=B/2*(R-r)./(r.*sind(phi));
F=2/pi*acos(exp(-f));

c=(8*pi*R.*F.*a.*lambda_r.*(sind(phi)).^2)./((1-a)*B*lambda.*Cn);
beta= (phi-alpha_d);

end