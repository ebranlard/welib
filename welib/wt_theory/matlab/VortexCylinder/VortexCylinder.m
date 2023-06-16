InitClear
%% for r=0
U0=10;
Uw=0.5*U0;
gamma=U0-Uw;
R=10;
D=2*R;
nD=2.5;
r=0;
x=linspace(-nD*D,nD*D,100);
ui=x*0;
for i=1:length(x)
    F = @(z0,psi)R*gamma/(4*pi)*(R-r*cos(psi))./((x(i)-z0).^2 + R^2+r^2-2*r*R*cos(psi)).^(3/2);
    ui(i) = dblquad(F,0,100*D,0,2*pi);
end
U=U0-ui;




figure
subplot(2,1,1)
hold on
ylim([0,2])
plot(x/D,ui/((U0-Uw)/2),'k','LineWidth',2)
plot([0 0],[0 2] ,'-','Color',[0.5,0.5,0.5],'LineWidth',2)
set(gca,'YTick',0:1:2)
[~,~] = format_ticks(gca,'',{'0','U_i','2U_i'})
%set(gca,'YTickLabel',{'0','Ui','2Ui'})
xlim([-2.5 2.5])
grid on 
box on
%xlabel('z/D [.]')
%ylabel('u_i(z) [m/s]')

% 
 subplot(2,1,2)
hold on
ylim([0,1])
xlim([-2.5 2.5])
plot(x/D,(U-Uw)/(U0-Uw),'k','LineWidth',2)
plot([0 0],[0 1] ,'-','Color',[0.5,0.5,0.5],'LineWidth',2)
set(gca,'YTick',0:0.5:1)
set(gca,'YTickLabel',{'U_w','U','U_0'});
grid on 
box on
[~,~] = format_ticks(gca,'',{'U_w','U','U_0'});
text(-0.1,-0.2,'z/D [.]')
%xlabel('')
%ylabel('U [m/s]')
title('VortexCylinderUUi')



%% For one r
U0=10;
Uw=1/3*U0;
gamma=U0-Uw;
R=40;
D=2*R;
nD=2.5;
r=0*R;
z=linspace(-nD*D,nD*D,20);
ui=z*0;


for i=1:length(z)
    F = @(z0,psi)R*gamma/(4*pi)*(R-r*cos(psi))./((z(i)-z0).^2 + R^2+r^2-2*r*R*cos(psi)).^(3/2);
    ui(i) = dblquad(F,0,100*D,0,2*pi);
end
U=U0-ui;


figure
hold on
ylim([0,1])
plot(z/D,U/U0,'k','LineWidth',2)
plot([min(z) max(z)]/D,[2/3 2/3] ,'k-.')
plot([min(z) max(z)]/D,[1/3 1/3] ,'k-.')
plot([0 0],[0 1] ,'-','Color',[0.5,0.5,0.5],'LineWidth',2)
grid on 
box on
xlabel('z/D [.]')
ylabel('U/U_0 [.]')
title('VortexCylinderU')

%% For lot's of r 3D
U0=10;
Uw=1/3*U0;
gamma=U0-Uw;
R=40;
D=2*R;
nD=2.5;
r=0*R;
z=linspace(-nD*D,nD*D,40);
r=linspace(0,2*R,30);
ui=zeros(length(r),length(z));

for k=1:length(r)
    for i=1:length(z)
        F = @(z0,psi)R*gamma/(4*pi)*(R-r(k)*cos(psi))./((z(i)-z0).^2 + R^2+r(k)^2-2*r(k)*R*cos(psi)).^(3/2);
        ui(k,i) = dblquad(F,0,100*D,0,2*pi);
    end
end
U=U0-ui;
%%
[X,Y]=meshgrid(z/D,r/R);
surf(X,Y,U/U0,'FaceColor','none')
xlabel('z/D [.]')
ylabel('r/R [.]')
zlabel('U/U_0 [.]')
xlim([-2.5 2.5])
zlim([-0 1.5])
ylim([0 2])
title('Vortex Cylinder 3D')
setFigureTitle(0)
%% for  3r in  2D
U0=10;
Uw=1/3*U0;
gamma=U0-Uw;
R=40;
D=2*R;
nD=2.5;
r=0*R;
z=linspace(-nD*D,nD*D,40);
r=[0 0.5 0.9]*R;
ui=zeros(length(r),length(z));

for k=1:length(r)
    for i=1:length(z)
        F = @(z0,psi)R*gamma/(4*pi)*(R-r(k)*cos(psi))./((z(i)-z0).^2 + R^2+r(k)^2-2*r(k)*R*cos(psi)).^(3/2);
        ui(k,i) = dblquad(F,0,100*D,0,2*pi);
    end
end
U=U0-ui;

figure
hold on
ylim([0,1])
plot(z/D,U(1,:)/U0,'k','LineWidth',2)
plot(z/D,U(2,:)/U0,'k--','LineWidth',2)
plot(z/D,U(3,:)/U0,'k-.','LineWidth',2)
plot([min(z) max(z)]/D,[2/3 2/3] ,'k-.')
plot([min(z) max(z)]/D,[1/3 1/3] ,'k-.')
plot([0 0],[0 1] ,'-','Color',[0.5,0.5,0.5],'LineWidth',2)
grid on 
box on
legend('r/R=0','r/R=0.5','r/R=0.9')
xlabel('z/D [.]')
ylabel('U/U_0 [.]')
title('VortexCylinderU')


%% At the rotor (z=0)
U0=10;
Uw=1/3*U0;
gamma=U0-Uw;
R=40;
D=2*R;
nD=2.5;
r=linspace(0,0.8*R,20);
z=0
ui=r*0;


for i=1:length(r)
    F = @(z0,psi)R*gamma/(4*pi)*(R-r(i)*cos(psi))./((z0).^2 + R^2+r(i)^2-2*r(i)*R*cos(psi)).^(3/2);
    ui(i) = dblquad(F,0,100*D,0,2*pi);
end
U=U0-ui;


figure
hold on
%ylim([0,1])
plot(r/R,U/U0,'k','LineWidth',2)
grid on 
box on
xlabel('r/R [.]')
ylabel('U/U_0 [.]')
title('VortexCylinderURotor')