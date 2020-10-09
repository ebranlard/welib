%% Initialization for the BEM code
InitClear;
InitDefault;
path(path,'./f_format')
path(path,'./f_aeroelastic')
path(path,'./f_optimal')

%
% Files={'data/test/NTK500_pitch_structure.htc','data/test/data/NewBlade.pc','data/test/data/NewBlade82.ae','data/NORDTANK_Spec.dat'};
% fInit('hawc',Files)

% 
% Files={'data/Tjaere/Tjaere_BladeGeometry.dat','data/Tjaere/Tjaere_BladeProfiles.dat','data/Tjaere/Tjaere_Spec.dat'};
% fInit('flex',Files)

% Files={'data/NORDTANK_BladeGeometry.dat','data/NORDTANK_BladeProfiles.dat','data/NORDTANK_Spec.dat'};
% fInit('flex',Files)

%  
% Files={'data/WTperf_aerodyn/Test03_CART3.wtp'};
% Files={'data/WTperf_aerodyn/Test02_AWT27.wtp'};
% fInit('wtperf',Files)




%%
R=Rotor.R;
r=Rotor.r;
% 
 lambda=8;
 alpha_d=6;  % alpha(17)=6 deg
 [a aprime phi c beta]=getOptimizedParameters(alpha_d,lambda,Rotor.nB);
% % 
 Rotor.chord=c;
 Rotor.beta=beta;

%% calling the BEM code 
Algo.nbIt=500;
Algo.alphaCrit=0.001;
Algo.relaxation=1;
Algo.correction='Spera';


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
%%

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

%ylim([5.98 6.1])

xlabel('r/R [.]')
ylabel('Angle of attack \alpha [deg]')
legend('BEM code','Optimum')


