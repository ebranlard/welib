%% Initialization for the BEM code
loading;
B=3;
R=3;
rhub=0.1;
lambda=8;
alpha_d=6;  % alpha(17)=6 deg
r= [0.3 0.5 1.0 1.5 2.0 2.5 2.6 2.7 2.8 2.9 2.95 ];
%size of each element dr
l=([rhub r]+[r R])/2;l(1)=rhub;l(length(l))=R;
dr=diff(l);     
x=r/R*lambda;   % this is actually lambda_r
[a aprime phi c beta ]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x);

%% calling the BEM code 
V0=8;
lambda=8;
Param.nbIt=5000;
Param.alphaCrit=0.05;
Param.relaxation=1;
Param.correction='none';
[BEM CP CT CQ] = BEMfunction(Param,Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho);

%% plotting
figure(21)
hold on
plot(r/R,lasts(BEM.a))
plot(r/R,a,'k+')
grid
box
xlabel('r/R [.]')
ylabel('Axial induction factor (a) [.]')
legend('BEM code','Optimum')

figure(22)
hold on
plot(r/R,lasts(BEM.phi)*180/pi)
plot(r/R,phi*180/pi,'k+')
grid
box
xlabel('r/R [.]')
ylabel('Flow angle (phi) [deg]')
legend('BEM code','Optimum')


figure(23)
hold on
plot(r/R,lasts(BEM.aprime))
plot(r/R,aprime,'k+')
grid
box
xlabel('r/R [.]')
ylabel('Tangential induction factor (a'') [.]')
legend('BEM code','Optimum')

%%

figure(24)
hold on
grid
box
plot(r/R,lasts(BEM.alpha))
plot(r/R,alpha_d*ones(1,11),'k+')
ylim([5.98 6.1])

xlabel('r/R [.]')
ylabel('Angle of attack \alpha [deg]')
legend('BEM code','Optimum')


