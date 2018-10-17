function [Vrel, FWind, TauWind, Ct ] = fWindCalc(t,q, Uw,zHub )

x=q(1); theta=q(2); xDot=q(3); thetaDot=q(4);

area = 0.25*pi*126^2;   % rotor area [m2]
Vrated = 11.4;          % [m/s]
rho = 1.29;             % air density [kg/m3]


Vrel = Uw - (xDot+thetaDot*zHub);

if Vrel <= Vrated
    Ct = 0.75;
else
    Ct = 0.75*exp(-0.25*(Vrel-Vrated)^.86);
end
if length(q)>4
    FWind = 0.5*rho*area*q(5)*Vrel*abs(Vrel);
else
    FWind = 0.5*rho*area*Ct*Vrel*abs(Vrel);
end
TauWind = FWind * zHub;

end

