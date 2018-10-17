function Vpoint=getPointIncomingWindOld(P, Tower, Angles, Model,VelocityParams)
yaw=Angles.yaw;
psi=Angles.psi;
Vhub=VelocityParams.V0;
nu=VelocityParams.nu;
Vpoint=[0 0 0];

V0=norm(Vhub); 


if isequal(Model,'Constant')
    Vpoint=Vhub;
end

if isequal(Model,'PowerLaw') || isequal(Model,'TowerEffectAndPowerLaw')
    Vpoint=Vhub*(P(1)/60).^nu;
    V0=norm(Vpoint);
end
if isequal(Model,'TowerEffect') || isequal(Model,'TowerEffectAndPowerLaw')
    if psi>80 & psi<280
          a=Tower.r1-(P(1)-49)*(Tower.r1/Tower.H1);
    elseif psi>120 & psi<240
          a=2.375-(P(1)-28)*7.8125e-3;
    else
          a=0;
    end
    % projection of 0P on the ground
    r=sqrt(P(3)^2+P(2)^2);
    cosT=P(3)/r;
    sinT=-P(2)/r;
    %Polar velocity
    Vr=V0*(1-(a/r)^2)*cosT;
    V_theta=-V0*(1+(a/r)^2)*sinT;
    %Cartesian Velocity
    % Wrong in case of slanted flow!!!!!!!!!1
    Vpoint(3)=Vr*cosT-V_theta*sinT;
    Vpoint(2)=-Vr*sinT-V_theta*cosT;
end

end


