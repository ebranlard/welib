function fPlotCompare(t1,q1,t2,q2,leg1,leg2,titl);
sty1='-';
sty2='-';
colr1=[0.35 0.35 0.35];
colr1=fColrs(3);
colr2=[0 0 0.7 ];
figure,hold all,box on,grid on
plot(t1, q1(:,1),sty1,'Color',colr1)
plot(t2, q2(:,1),sty2,'Color',colr2,'LineWidth',1)
ylabel('Surge x0 [m]')
xlabel('t [s]')
legend(leg1,leg2)
title([titl 'PlotCompSurge']);


figure,hold all,box on,grid on
plot(t1, q1(:,2)*180/pi,sty1,'Color',colr1)
plot(t2, q2(:,2)*180/pi,sty2,'Color',colr2,'LineWidth',1)
ylabel('Pitch \theta [deg]')
xlabel('t [s]')
legend(leg1,leg2)
title([titl 'PlotCompPitch']);

