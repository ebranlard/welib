%% Checking theory
r=linspace(0,3,20);
beta=0+zeros(1,length(r));
EI1=2400;
EI2=24000;
V=20;
nu=0;   
n=length(r);
%Assuming constant load distribution
py=ones(1,n)*100;
pz=ones(1,n)*100;

%getting uz and uy
[uy uz]=fDeflectionFromLoad(py,pz,beta,nu,r,EI1,EI2);

R=3;
xi=(R-r)/R;
%uanaly=py*(R^4)./(24*EI2).*(xi.^4-4.*xi+3);
uanalz=pz*(R^4)./(24*EI1).*(xi.^4-4.*xi+3);;

uanaly=py./(24*EI2).*r.^2.*(2*R*(3*R-2*r)+r.^2);

figure(1)
hold on
plot(r,uanaly,'k+')
plot(r,uy,'b')
grid 
box
legend('Beam theory', 'Algorithm','Location','NorthWest')
xlabel('r [m]')
ylabel('uy [m]')

figure(2)
hold on
plot(r,uanalz,'k+')
plot(r,uz,'b')
grid
box
legend('Beam theory', 'Algorithm','Location','NorthWest')
xlabel('r [m]')
ylabel('uz [m]')
