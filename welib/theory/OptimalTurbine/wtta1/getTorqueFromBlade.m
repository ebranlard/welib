function Q=getTorqueFromBlade(r0,Pt0,R)
n=length(r0);
r=zeros(1,n+1);
Pt=zeros(1,n+1);
r(1:n)=r0;
r(end)=R;
Pt(1:n)=Pt0;

Q=trapz(r,r.*Pt);

Mtot=0;
for i=1:length(r0)
    Mtot=Mtot+ 1/3*(Pt(i+1)-Pt(i))/(r(i+1)-r(i))  *(r(i+1)^3-r(i)^3)  + 1/2* (Pt(i)*r(i+1) -Pt(i+1)*r(i))/(r(i+1)-r(i)) *(r(i+1)^2-r(i)^2);
end
%Q=Mtot;

end