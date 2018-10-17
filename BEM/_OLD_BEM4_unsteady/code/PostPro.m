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
% 
% figure(31)
% hold all
% plot(t,x(1,:))
% plot(t,v(1,:))
% plot(t,a(1,:))
% xlabel('Time[s]')
% title('M_1')
% grid on
% legend('Displacement [m]','Velocity [m/s]', 'Acceleration [m/s^2]','location','north')
% xlim([0 T0])

% 
% 
% %% plots in term of theta
% theta=x(2,:)*180/pi;
% Theta0=180;
% imax=getMax(theta);
% imax=imax(1)+1
% It=1:imax;
% It2=imax:(2*imax);
% theta2=theta((It2));
% theta=theta((It));
% 
% close all
% setFigureWidth('18')
% setFigureHeight('11')
% figure(30)
% hold all
% plot(theta,x(2,It))
% plot(theta,v(2,It))
% plot(theta,a(2,It))
% plot(theta2,v(2,It2),'-g')
% xlabel('\theta [^o]')
% title('Rodt')
% grid on
% legend('Displacement [rad]','Velocity [rad/s]', 'Acceleration [rad/s^2]','location','northEast')
% xlim([0 Theta0])
% 
% figure(31)
% hold all
% plot(theta,x(1,It))
% plot(theta,v(1,It))
% plot(theta,a(1,It))
% plot(theta2,v(1,It2),'-g')
% xlabel('\theta [^o]')
% title('M_1t')
% grid on
% legend('Displacement [m]','Velocity [m/s]', 'Acceleration [m/s^2]','location','northEast')
% xlim([0 Theta0])
% 
% 

