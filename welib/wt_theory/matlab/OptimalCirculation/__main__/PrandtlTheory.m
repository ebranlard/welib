InitClear;
setFigureWidth('0');

theta=linspace(0,15*pi,10000);
t=cos(theta)+1i*sin(theta);

a=5;

jouk=0.5*(t+1./t);
joukre=real(jouk);
% joukre=jouk;

z=a/pi*log(joukre);
%
figure(12)
clf
subplot(1,3,1)
plot(t)
subplot(1,3,2)
plot(jouk)
ylim([-1 1])
subplot(1,3,3)
plot(z)
ylim([-12 12])


%%
figure
hold all
n=15;
epsilon=0.5;
vr=linspace(1+epsilon,5,n)

X=[];
Y=[];
U=[];
V=[];
for i=1:n;
r=vr(i);
theta=linspace(-pi,pi,75);
Z=r*exp(1i*theta);


W=-1i*(Z.^2+1)./(Z.^2-1);
z=a/pi*log(0.5*(Z+1./Z));
x=real(z);
y=imag(z);

u=real(W);
v=-imag(W);
%
X=[X x];
Y=[Y y];
U=[U u];
V=[V v];

% quiver(-y,x,-v,u,'Color',fColrs(i))
% xlabel('y')
% ylabel('x')
end
plot([-1 -1],[-100 0],'k','LineWidth',4.2)
plot([1 1],[-100 0],'k','LineWidth',4.2)
plot([0 0],[-100 0],'k','LineWidth',4.2)
quiver(Y/a,X/a,V,U,'Color','k')
xlabel('$y/a$ [.]')
ylabel('$x/a$ [.]')
xlim([-1 1])
box on
ylim([-0.4 0.4])
title('PrandtlVelocityField')



%%
x=linspace(-1000,0,15000);
Phi=a/pi*acos(exp(x*pi/a))/a;
Phi2=-a/pi*acos(exp(x*pi/a))/a;


figure
hold all
plot(Phi,x/a,'k')
plot(Phi2,x/a,'k')
plot(Phi-1,x/a,'k')
plot(Phi2-1,x/a,'k')
plot(Phi+1,x/a,'k')
plot(Phi2+1,x/a,'k')
plot(Phi-2,x/a,'k')
plot(Phi2+2,x/a,'k')
xlim([-2 2])
ylim([-1.5 0])
xlabel('$\Phi/a$ [m/s]')
ylabel('$x/a$ [.]')
box on
title('PrandtlPotential')
setFigureLatex(1)

%%
trapz(x,Phi)/2
%%
figure
hold all
n=15;
vk=linspace(-0.5*a,a,n);

X=[];
Y=[];
U=[];
V=[];
for i=1:n;
x=ones(1,150)*vk(i);
y=linspace(-3*a,3*a,150);
z=x+1i*y    ;
Z=-exp(pi*z/(a))-1i*sqrt(1-exp(pi*z/(a)).^2);


W=-1i*(Z.^2+1)./(Z.^2-1);


u=real(W);
v=-imag(W);
%
X=[X x];
Y=[Y y];
U=[U u];
V=[V v];

% quiver(-y,x,-v,u,'Color',fColrs(i))
% xlabel('y')
% ylabel('x')
end
plot([-1 -1],[-100 0],'k','LineWidth',4.2)
plot([1 1],[-100 0],'k','LineWidth',4.2)
plot([0 0],[-100 0],'k','LineWidth',4.2)
quiver(Y/a,X/a,V,U,'Color','k')
xlabel('$y/a$ [.]')
ylabel('$x/a$ [.]')
xlim([-1 1])
box on
ylim([-0.4 0.4])
title('PrandtlVelocityField')


%%
InitClear
Vr=linspace(0,1,200);
VB=[20 3 2 1];
Vsty={'-','--','-.',':'};
figure(1),box on
set(1,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|--|-.|:');
legs=cell(1,length(VB));
hold all
for iB=1:length(VB)
    B=VB(iB);
    F=fTipLossPrandtl(7,1,B,Vr,1);
    plot(Vr,F)
    legs{iB}=sprintf('$B = %d$',B);
end
grid on
xlabel('r/R [.]')
ylabel('F [.]')
title('TiplossPrandtlB')
legend(legs,'Location','SouthWest')
setFigureLatex(1);
xlim([0.3 1])


