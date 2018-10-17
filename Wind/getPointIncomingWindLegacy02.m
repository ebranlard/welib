function Vpoint=getPointIncomingWindLegacy02(P,psi,WT,Wind,Algo)
Model=Wind.Model;
Vhub=Wind.V0;
nu=Wind.nu;
Vpoint=[0;0;0];
V0=norm(Vhub); 


if length(findstr('Constant',Model))>0
    Vpoint=Vhub;
end

if length(findstr('PowerLaw',Model))>0
    Vpoint=Vhub*(P(1)/60).^nu;
    V0=norm(Vpoint);
end
if length(findstr('TowerEffect',Model))>0
    if psi>80 & psi<280
          a=WT.Tower.r1-(P(1)-49)*(WT.Tower.r1/WT.Tower.H1);
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

if length(findstr('Stochastic',Model))>0
    Vpoint=Vpoint+Aero.Wind.fTurb(P(3),-P(2),P(1),Algo.t);
end


end


