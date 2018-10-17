function Q=getTorqueFromBlade(r,Pt)
R=[r ;30.56];
Pt=[Pt.*r;0];
Q=trapz(R,Pt);
end