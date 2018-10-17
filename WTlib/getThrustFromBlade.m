function T=getThrustFromBlade(r0,Pn0,R)
n=length(r0);
r=zeros(1,n+1);
Pn=zeros(1,n+1);
r(1:n)=r0;
r(end)=R;
Pn(1:n)=Pn0;

T=trapz(r,Pn);
end