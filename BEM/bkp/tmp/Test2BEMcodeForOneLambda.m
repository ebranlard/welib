%% Initialization for the BEM code
InitClear
InitRotor
InitSystemProperties
R=Rotor.R;
r=Rotor.r;

lambda=8;
alpha_d=6;  % alpha(17)=6 deg
[a aprime phi c beta]=getOptimizedParameters(alpha_d,lambda,Rotor.nB);

Rotor.chord=c;
Rotor.beta=beta;

%% calling the BEM code 
BEMParam.nbIt=500;
BEMParam.alphaCrit=0.01;
BEMParam.relaxation=1;
BEMParam.correction='none';


V0=8;
Aero.Wind.V0=[0; 0; V0];
Rotor.Omega=lambda/R*V0;
% calling the BEM code for a given lambda and wind speed
BEM = fBEMsteady();
%% plotting
figure(21)
hold on
plot(r/R,BEM.a)
plot(r/R,a,'k+')
grid
box
xlabel('r/R [.]')
ylabel('Axial induction factor (a) [.]')
legend('BEM code','Optimum')

figure(22)
hold on
plot(r/R,BEM.phi)
plot(r/R,phi,'k+')
grid
box
xlabel('r/R [.]')
ylabel('Flow angle (phi) [deg]')
legend('BEM code','Optimum')


figure(23)
hold on
plot(r/R,BEM.aprime)
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
plot(r/R,BEM.alpha)
plot(r/R,alpha_d*ones(1,Rotor.ne),'k+')
ylim([5.98 6.1])

xlabel('r/R [.]')
ylabel('Angle of attack \alpha [deg]')
legend('BEM code','Optimum')


