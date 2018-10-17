function Q=getTorqueFromBlade(r,Pt)
global Rotor
%R=[r ;Rotor.R];
R=[r ;30.56];
Pt=[Pt;0];
% trapz is not accurate due to the curvature
%  Q=trapz(R,Pt.*R);
Q=0;
for j=1:(length(r))
    A(j)=(Pt(j+1)-Pt(j))/(R(j+1)-R(j));
    B(j)=(((Pt(j)*R(j+1))-(Pt(j+1))*R(j)))/(R(j+1)-R(j));
    dQ(j)=(A(j)/3)*((R(j+1))^3-(R(j))^3)+0.5*(B(j))*((R(j+1))^2-(R(j))^2);
    Q = Q + dQ(j);
end


end