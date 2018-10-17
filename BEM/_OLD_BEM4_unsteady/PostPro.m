setFigurePath('../figs/');
setFigureTitle(0);
%% plots in terms of t
close all
%setFigureWidth('18')
%setFigureHeight('11')
%%


%%
setFigureWidth('24')
setFigureHeight('20')
figure
subplot(4,1,1);plot(t,Thrust/1000);grid on;xlabel('t [s]');ylabel('Thrust [kN]');xlim(XLIM);
subplot(4,1,2);plot(t,Torque/1000);grid on;xlabel('t [s]');ylabel('Torque [kNm]');xlim(XLIM);
subplot(4,1,3);plot(t,MG/1000);grid on;xlabel('t [s]');ylabel('MG [kNm]');xlim(XLIM);
subplot(4,1,4);plot(t,Power/1000);grid on;xlabel('t [s]');ylabel('Power [kW]');xlim(XLIM);
title(strcat('Aero',flag));

%%
setFigureWidth('24')
setFigureHeight('15')
figure
hold all
subplot(3,1,1);plot(t,x(1,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Position x [m]');
subplot(3,1,2);plot(t,v(1,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed dx/dt [m/s]');
subplot(3,1,3);plot(t,a(1,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Acc.  d^2x/dt^2 [m/s^2]');
title(strcat('Nacelle',flag))
%legend('Displacement [rad]','Velocity [rad/s]', 'Acceleration [rad/s^2]','location','north')
%%
setFigureWidth('24')
setFigureHeight('15')
figure;hold all
subplot(3,1,1);plot(t,sin(x(2,:)));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Position sin(\theta) [rad]');
subplot(3,1,2);plot(t,v(2,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed d\theta/dt [rad/s]');
subplot(3,1,3);plot(t,a(2,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Acc.  d^2\theta/dt^2 [rad/s^2]');
title(strcat('Rotor',flag))
%legend('Displacement [rad]','Velocity [rad/s]', 'Acceleration [rad/s^2]','location','north')
%%
setFigureWidth('24')
setFigureHeight('15')
figure;hold all
subplot(3,1,1);plot(t,sin(x(3,:)));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Position sin(\nu) [rad]');
subplot(3,1,2);plot(t,v(3,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed d\nu/dt [rad/s]');
subplot(3,1,3);plot(t,a(3,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Acc.  d^2\nu/dt^2 [rad/s^2]');
title(strcat('Generator',flag))
grid on
%%

% 
b=Rotor.Blade;
r=Rotor.r;
uz=zeros(3,nt,ne);
uy=zeros(3,nt,ne);
for idB=1:3
    for i=1:nt
        uz(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,3)+x(idB*3+2,i)*b.eigen1e(:,3)+x(idB*3+3,i)*b.eigen2f(:,3);
        uy(idB,i,:)=x(idB*3+1,i)*b.eigen1f(:,2)+x(idB*3+2,i)*b.eigen1e(:,2)+x(idB*3+3,i)*b.eigen2f(:,2);
    end
end
    %%

% % 
% for i=1:nt
%     if(i>1 && t(i)==0) break;
%     end
%     figure(5)
%     clf
%     hold all
%     grid on
%     plot3(r,squeeze(uy(1,i,:)),squeeze(uz(1,i,:)))
%     %    plot3(r,squeeze(uy(1,i,:)),squeeze(uz(1,i,:)))    
%     xlim([6 31])
%     ylim([min(min(min(uy))) max(max(max(uy)))])
%     zlim([min(min(min(uz))) max(max(max(uz)))])
%     view(3)
%     pause(0.01)
%     
% end    


%%
I=floor(10/dt):length(t);
% I=1:floor(10/dt);
% I=1:length(t);

figure
hold all
plot(mod(x(2,I)*180/pi,360),uy(1,I,end),'.')
plot(mod(x(2,I)*180/pi,360),uy(2,I,end),'.')
plot(mod(x(2,I)*180/pi,360),uy(3,I,end),'.')
title('Uy')
%%
figure
hold all
plot(x(2,I)*180/pi-x(2,1)*180/pi,uy(1,I,end),'-')
plot(x(2,I)*180/pi-x(2,1)*180/pi,uy(2,I,end),'-')
plot(x(2,I)*180/pi-x(2,1)*180/pi,uy(3,I,end),'-')
title('Uy')

%%
figure
hold all
plot(mod(x(2,I)*180/pi,360),uz(1,I,end),'.')
plot(mod(x(2,I)*180/pi,360),uz(2,I,end),'.')
plot(mod(x(2,I)*180/pi,360),uz(3,I,end),'.')
title('Uz')

%%

[x1,y1]=bin(mod(x(2,I)*180/pi,360),uy(1,I,end),30);
[x2,y2]=bin(mod(x(2,I)*180/pi+120,360),uy(2,I,end),30);
[x3,y3]=bin(mod(x(2,I)*180/pi+240,360),uy(3,I,end),30);
figure
plot(x1,y1,x2,y2,x3,y3)
title('Uy')
%%
[x1,y1]=bin(mod(x(2,I)*180/pi,360),uz(1,I,end),30);
[x2,y2]=bin(mod(x(2,I)*180/pi+120,360),uz(2,I,end),30);
[x3,y3]=bin(mod(x(2,I)*180/pi+240,360),uz(3,I,end),30);
figure
plot(x1,y1,x2,y2,x3,y3)
title('Uz')
%%
figure
hold all
plot(t(I),uz(1,I,end))
plot(t(I),uz(2,I,end)+1)
plot(t(I),uz(3,I,end)+2)
title('Flap motion')
%%
figure
hold all
plot(t(I),Flap(I,1))
plot(t(I),Flap(I,2))
plot(t(I),Flap(I,3))
title('Flapwise moment')
%%
figure
hold all
plot(t,Edge(:,1))
plot(t,Edge(:,2))
plot(t,Edge(:,3))
title('Edgewise moment')


%%
setFigureWidth('24')
setFigureHeight('15')
figure;hold all
subplot(3,1,1);plot(t,sin(x(4,:)));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Position [arb]');
subplot(3,1,2);plot(t,v(4,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed [arb]');
subplot(3,1,3);plot(t,a(4,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Acc [arb]');
title(strcat('Blade11flap',flag))
grid on

%%
setFigureWidth('24')
setFigureHeight('15')
figure;hold all
subplot(3,1,1);plot(t,sin(x(5,:)));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Position [arb]');
subplot(3,1,2);plot(t,v(5,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Speed [arb]');
subplot(3,1,3);plot(t,a(5,:));xlim(XLIM);grid on;xlabel('t [s]');ylabel('Acc. [arb]');
title(strcat('Blade11edge',flag))
grid on
