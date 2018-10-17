%% Systems initialization
InitClear;
InitSystemProperties;  % Aero, Rotor, Nacelle, Shaft, Generator, Tower are definied
flag='12';

%% Main Parameters
global Algo
dt=0.02; % time step [s]
tmax=30; % simulation time [s]

% Algorithm switches
Algo.YawModel=0;
Algo.dynastall=1;
Algo.dt=dt;
Algo.damp=1*10^6;
Algo.TwoDOF=0;
Algo.Weight=0;
Algo.DontRotate=0;
Aero.Wind.fV0=@(t)[0; 0;10 ] ; % wind speed at hub height
% Controller.fpitch=@(t) 0;

% Initial condition
a0=[ -739.8815e-006;    23.0133e-003;  -104.3992e-003; zeros(9,1)]*1;
v0=[ -2.9467e-003;  2.2997e+000;  30.4165e-003;zeros(9,1)]*1;
x0=[ 78.4824e-003; 36.8422e+000; -8.3656e-003; zeros(9,1)]*1;


%% System resolution with Runge Kutta nystrom scheme
Algo.Weight=0;
Aero.Wind.Model='PowerLawConstant';
Runge

%% 
b=Rotor.Blade;r=Rotor.r;uz=zeros(3,nt,ne);uy=zeros(3,nt,ne);
for idB=1:3
    for i=1:nt
        uz(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,3)+x(idB*3+2,i)*b.eigen1e(:,3)+x(idB*3+3,i)*b.eigen2f(:,3);
        uy(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,2)+x(idB*3+2,i)*b.eigen1e(:,2)+x(idB*3+3,i)*b.eigen2f(:,2);
    end
end
xpl=x;
vpl=v;
Flappl=Flap;
Edgepl=Edge;
uypl=uy;
uzpl=uz;
%%
setFigureWidth('0')
setFigureHeight('0')
I=floor(10/dt):length(t);
I=913:length(t);
% I=1:floor(10/dt);
% I=1:length(t);
[x1,y1]=bin(mod(xpl(2,I)*180/pi,360),uypl(1,I,end),30);
[x2,y2]=bin(mod(xpl(2,I)*180/pi,360),uypl(2,I,end),30);
[x3,y3]=bin(mod(xpl(2,I)*180/pi,360),uypl(3,I,end),30);

figure
hold on
plot(mod(xpl(2,I)*180/pi,360),uypl(1,I,end),'b.')
plot(mod(xpl(2,I)*180/pi,360),uypl(2,I,end),'g.')
plot(mod(xpl(2,I)*180/pi,360),uypl(3,I,end),'r.')
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k','LineWidth',1)
% ylim([0 0.25])
xlim([0 360])
grid on
title('UyPowerLaw')
ylabel('Tip Deflection U_y [m]')
xlabel('Azimuthal position \psi [deg]')
%%
[x1,y1]=bin(mod(xpl(2,I)*180/pi,360),uzpl(1,I,end),30);
[x2,y2]=bin(mod(xpl(2,I)*180/pi,360),uzpl(2,I,end),30);
[x3,y3]=bin(mod(xpl(2,I)*180/pi,360),uzpl(3,I,end),30);

figure
hold all
plot(mod(xpl(2,I)*180/pi,360),uzpl(1,I,end),'b.')
plot(mod(xpl(2,I)*180/pi,360),uzpl(2,I,end),'g.')
plot(mod(xpl(2,I)*180/pi,360),uzpl(3,I,end),'r.')
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k','LineWidth',1)
title('UzPowerLaw')
ylabel('Tip Deflection U_z [m]')
xlabel('Azimuthal position \psi [deg]')

xlim([0 360])
% ylim([1.05 1.4])
grid on

%% PostProcessing
Algo.Weight=1; 
Aero.Wind.Model='Constant';
Runge

%% 

b=Rotor.Blade;r=Rotor.r;uz=zeros(3,nt,ne);uy=zeros(3,nt,ne);
for idB=1:3
    for i=1:nt
        uz(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,3)+x(idB*3+2,i)*b.eigen1e(:,3)+x(idB*3+3,i)*b.eigen2f(:,3);
        uy(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,2)+x(idB*3+2,i)*b.eigen1e(:,2)+x(idB*3+3,i)*b.eigen2f(:,2);
    end
end
xw=x;
vw=v;
Flapw=Flap;
Edgew=Edge;
uyw=uy;
uzw=uz;
%%
I=913:length(t);
[x1,y1]=bin(mod(xw(2,I)*180/pi,360),uyw(1,I,end),30);
[x2,y2]=bin(mod(xw(2,I)*180/pi,360),uyw(2,I,end),30);
[x3,y3]=bin(mod(xw(2,I)*180/pi,360),uyw(3,I,end),30);

figure
hold on
plot(mod(xw(2,I)*180/pi,360),uyw(1,I,end),'b.')
plot(mod(xw(2,I)*180/pi,360),uyw(2,I,end),'g.')
plot(mod(xw(2,I)*180/pi,360),uyw(3,I,end),'r.')
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k','LineWidth',1)
title('UyWeight')
ylim([0 0.25])
legend()
xlim([0 360])
grid on
ylabel('Tip Deflection U_y [m]')
xlabel('Azimuthal position \psi [deg]')

%%
[x1,y1]=bin(mod(xw(2,I)*180/pi,360),uzw(1,I,end),30);
[x2,y2]=bin(mod(xw(2,I)*180/pi,360),uzw(2,I,end),30);
[x3,y3]=bin(mod(xw(2,I)*180/pi,360),uzw(3,I,end),30);

figure
hold on
plot(mod(xw(2,I)*180/pi,360),uzw(1,I,end),'b.')
plot(mod(xw(2,I)*180/pi,360),uzw(2,I,end),'g.')
plot(mod(xw(2,I)*180/pi,360),uzw(3,I,end),'r.')
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k','LineWidth',1)
xlim([0 360])
ylim([1.05 1.4])
grid on
title('UzWeight')
ylabel('Tip Deflection U_z [m]')
xlabel('Azimuthal position \psi [deg]')




%%

%%

%%

%%

%%

%%

%%

%%

%%

%%

%%


%% PostProcessing
Algo.Weight=0; 
Aero.Wind.Model='ConstantTowerEffect';
Runge

%% 

b=Rotor.Blade;r=Rotor.r;uz=zeros(3,nt,ne);uy=zeros(3,nt,ne);
for idB=1:3
    for i=1:nt
        uz(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,3)+x(idB*3+2,i)*b.eigen1e(:,3)+x(idB*3+3,i)*b.eigen2f(:,3);
        uy(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,2)+x(idB*3+2,i)*b.eigen1e(:,2)+x(idB*3+3,i)*b.eigen2f(:,2);
    end
end
xw=x;
vw=v;
Flapw=Flap;
Edgew=Edge;
uyw=uy;
uzw=uz;
%%
I=floor(10/dt):length(t);
I=913:length(t);

[x1,y1]=bin(mod(xw(2,I)*180/pi,360),uyw(1,I,end),30);
[x2,y2]=bin(mod(xw(2,I)*180/pi,360),uyw(2,I,end),30);
[x3,y3]=bin(mod(xw(2,I)*180/pi,360),uyw(3,I,end),30);

figure
hold on
plot(mod(xw(2,I)*180/pi,360),uyw(1,I,end),'b.')
plot(mod(xw(2,I)*180/pi,360),uyw(2,I,end),'g.')
plot(mod(xw(2,I)*180/pi,360),uyw(3,I,end),'r.')
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k','LineWidth',1)
title('UyTower')
ylim([0 0.25])
legend()
xlim([0 360])
grid on
ylabel('Tip Deflection U_y [m]')
xlabel('Azimuthal position \psi [deg]')

%%
[x1,y1]=bin(mod(xw(2,I)*180/pi,360),uzw(1,I,end),30);
[x2,y2]=bin(mod(xw(2,I)*180/pi,360),uzw(2,I,end),30);
[x3,y3]=bin(mod(xw(2,I)*180/pi,360),uzw(3,I,end),30);

figure
hold on
plot(mod(xw(2,I)*180/pi,360),uzw(1,I,end),'b.')
plot(mod(xw(2,I)*180/pi,360),uzw(2,I,end),'g.')
plot(mod(xw(2,I)*180/pi,360),uzw(3,I,end),'r.')
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k','LineWidth',1)
xlim([0 360])
ylim([1.05 1.4])
grid on
title('UzTower')
ylabel('Tip Deflection U_z [m]')
xlabel('Azimuthal position \psi [deg]')



