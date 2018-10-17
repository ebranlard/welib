function T=getThrustFromBlade(r,Pn)
R=[r' 30.56];
Pn=[Pn' 0];
T=trapz(R,Pn);
end