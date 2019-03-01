function RES = fakeBEMClCd(Cl,Cd,r,c,a,aprime,B,R,lambda,V0,rho,phi0)
    global Rotor
    %variable needed  
    ne=length(r);
    omega=lambda*V0/R;
    sigma=c*B./(2*pi*r);
    %% parameters for each element

    RES.Cn=zeros(ne,1);
    RES.Ct=zeros(ne,1);
    RES.F=zeros(ne,1);
    %% loads paramters for each element
    RES.Vr=zeros(ne,1);
    RES.Ct=zeros(ne,1);
    RES.Cp=zeros(ne,1);
    RES.L=zeros(ne,1);
    RES.D=zeros(ne,1);
    RES.T=zeros(ne,1);
    RES.Q=zeros(ne,1);
    RES.P=zeros(ne,1);
   
    RES.Pn=zeros(ne,1);
    RES.Pt=zeros(ne,1);
    
     %step 2
    RES.a=a';
    RES.aprime=aprime';
     %step 2
    RES.phi=phi0;
    %step 3
    for e=1:ne
        %step 4
        RES.Cl = Cl;
        RES.Cd = Cd;
        %step 5
        f=B/2*(R-r(e))/(r(e)*sind(RES.phi(e)));
        RES.F(e)=2/pi*acos(exp(-f));

        %% loads
        Vr=V0*(1-RES.a(e))/sind(RES.phi(e));
        phi=RES.phi(e);
        %%% Step 5 : Aerodynamic forces PER LENGTH
        RES.L=0.5*rho*norm(Vr)^2*c(e)*RES.Cl;
        RES.D=0.5*rho*norm(Vr)^2*c(e)*RES.Cd;
        RES.Pn(e) = RES.L*cosd(phi) + RES.D*sind(phi);  %load normal to the rotor plane
        RES.Pt(e) = RES.L*sind(phi) - RES.D*cosd(phi);   %load tangential to the rotor plane
        RES.Cn(e)=RES.Cl*cosd(phi)+RES.Cd*sind(phi);
        RES.Ct(e)=RES.Cl*sind(phi)-RES.Cd*cosd(phi);        
    end %end for loop on elements
    BladeTorque=getTorqueFromBlade(r,RES.Pt,R);
    BladeThrust=getThrustFromBlade(r,RES.Pn,R);             
    Torque=B*BladeTorque;
    Power=Torque*omega;        
    RES.CP=Power/(0.5*rho*V0^3*pi*R^2);
    RES.CT=BladeThrust/(0.5*rho*V0^2*pi*R^2);
    RES.CQ=Torque/(0.5*rho*V0^2*pi*R^3);
end

