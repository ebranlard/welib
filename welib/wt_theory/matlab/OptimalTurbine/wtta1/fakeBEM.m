function [BEM CP CT CQ] = fakeBEM(Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho,phi)
    %variable needed  
    ne=length(r);
    omega=lambda*V0/R;
    sigma=c*B./(2*pi*r);
    %% matrices that will store the paramters value at each iterational step
    BEM.a=zeros(ne,1);
    BEM.aprime=zeros(ne,1);
    BEM.alpha=zeros(ne,1);
    BEM.phi=zeros(ne,1);
    BEM.Cl=zeros(ne,1);
    BEM.Cd=zeros(ne,1);
    BEM.Cn=zeros(ne,1);
    BEM.Ct=zeros(ne,1);
    BEM.F=zeros(ne,1);
    %% vectors that will store the loads paramters for each elements
    BEM.Vr=zeros(ne,1);
    BEM.Ct=zeros(ne,1);
    BEM.Cp=zeros(ne,1);
    BEM.L=zeros(ne,1);
    BEM.D=zeros(ne,1);
    BEM.T=zeros(ne,1);
    BEM.Q=zeros(ne,1);
    BEM.P=zeros(ne,1);
    for e=1:ne
        %initialization
        BEM.a(e,1)=a(e);
        BEM.aprime(e,1)=aprime(e);
        %step 2
        BEM.phi(e,1)=phi(e);
        %step 3
        BEM.alpha(e,1)=(BEM.phi(e,1)*180/pi)-(beta(e)); %deg
        %step 4
        id=whichvalue(Profile.alpha,BEM.alpha(e,1));
        BEM.Cl(e,1)=Profile.Cl(id);
        BEM.Cd(e,1)=Profile.Cd(id);
        %step 5
        f=B/2*(R-r(e))/(r(e)*sin(BEM.phi(e,1)));
        BEM.F(e,1)=2/pi*acos(exp(-f));
        BEM.Cn(e,1)=BEM.Cl(e,1)*cos(BEM.phi(e,1))+BEM.Cd(e,1)*sin(BEM.phi(e,1));
        BEM.Ct(e,1)=BEM.Cl(e,1)*sin(BEM.phi(e,1))-BEM.Cd(e,1)*cos(BEM.phi(e,1));
        %% loads
        BEM.Vr(e)=V0*(1-BEM.a(e,1))/sin(BEM.phi(e,1));
        BEM.L(e)=BEM.Cl(e,1)*0.5*rho*BEM.Vr(e)^2*c(e)*dr(e);
        BEM.D(e)=BEM.Cd(e,1)*0.5*rho*BEM.Vr(e)^2*c(e)*dr(e);
        BEM.T(e)=B*BEM.Cn(e,1)*0.5*rho*BEM.Vr(e)^2*c(e)*dr(e);
        BEM.Q(e)=B*BEM.Ct(e,1)*0.5*rho*BEM.Vr(e)^2*c(e)*dr(e)*r(e);
        BEM.P(e)=BEM.Q(e)*omega;
        BEM.Cp(e)=BEM.P(e)/(0.5*rho*V0^3*pi*R^2);
    end %end for loop on elements
    CP=sum(noNA(BEM.Cp));
    CT=sum(BEM.T)/(0.5*rho*V0^2*pi*R^2);
    CQ=sum(BEM.Q)/(0.5*rho*V0^2*pi*R^3);
end

