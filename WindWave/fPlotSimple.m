function fPlotSimple(t1,q1,titl);
sty1='-';
sty2='-';
colr1=[0.35 0.35 0.35];
colr1=fColrs(3);
colr2=[0 0 0.7 ];
figure,hold all,box on,grid on
plot(t1, q1(:,1),sty1,'Color',colr2)
ylabel('Surge x0 [m]')
xlabel('t [s]')
title([titl 'PloSimpSurge']);


figure,hold all,box on,grid on
plot(t1, q1(:,2)*180/pi,sty1,'Color',colr2)
ylabel('Pitch \theta [deg]')
xlabel('t [s]')
title([titl 'PlotSimpPitch']);

