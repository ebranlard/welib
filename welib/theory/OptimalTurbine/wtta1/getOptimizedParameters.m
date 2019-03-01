function [a aprime phi c beta]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x)
    Cd_opt=Profile.Cd(whichvalue(Profile.alpha,alpha_d));
    Cl_opt=Profile.Cl(whichvalue(Profile.alpha,alpha_d));
    a=zeros(1,length(x));
    for i=1:length(x)
        a(i)=fzero(@(aa) getRelationAX(aa,x(i)),0.3);
    end    
    aprime=(1-3*a)./(4*a-1);
    phi=atan( (1-a)./((1+aprime).*x) );  %rad
    Cn=Cl_opt*cos(phi)+Cd_opt*sin(phi);
    f=B/2*(R-r)./(r.*sin(phi));
    F=2/pi*acos(exp(-f));
    c=(8*pi*R.*F.*a.*x.*(sin(phi)).^2)./((1-a)*B*lambda.*Cn);
    beta= ((phi*180/pi)-alpha_d);
end