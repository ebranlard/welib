%%  Init
InitClear;
PATH.EXPANSION='/work/lib/WTTheory/WakeExpansion/';
require('EXPANSION','');

require('THEODORSEN','v02');
require('OPTIMCIRC','v01_nomex');


%% Params
CT=0.737;
lambda=7.06;
nB=3;

U0=10;% should not matter
R=1;  % should not matter

%% Different models
vz_barb=linspace(0,30,100); %z/R !!! not z/D

[ exp_franksen ] = fExpansionFranksen(CT,0.7,2,vz_barb);
[ exp_rathmann ] = fExpansionRathmann(CT,vz_barb);
[ exp_vortexrings ] = fExpansionVortexRings(CT,vz_barb); 

%% Theodorsen - Not iterative
[ exp_theodorsen_notiterative vz_bar] = fTheodorsenExpansionWrap(lambda,nB,U0,R,1/2*(1-sqrt(1-CT)),CT); % !!!! mean a, that's not really good

%% Theodorsen iterative
[w_bar l_bar Rw exp_theodorsen_iterative]=fTheodorsenFarWakeParams(CT,lambda,nB);

%% Theodorsen Using Goldstein's(Okulov) far wake param... Needs lower level interface fTheodorsenExpansion
[w_bar_gold l_bar_gold CTout G ]=fGoldsteinFarWakeParams(CT,lambda,nB);
Algo.nr=100;
Algo.ntheta=3000;
Algo.z_inf=100;
Algo.bWT=1;
vx=linspace(0,1,Algo.nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
[ ~, ExpFactor] = fTheodorsenExpansion(l_bar_gold,vx,G,nB,Algo.z_inf,Algo.ntheta,Algo.bWT );
exp_theodorsen_goldstein=1+ExpFactor*CT;


%% Plot
k=0.05;m=1;
figure,hold all,grid on,box on
plot(vz_barb/2,exp_franksen,'--')
plot(vz_barb/2,exp_rathmann,'--')
plot(vz_barb/2,exp_vortexrings,'k--')
% plot(vz_barb/2,k*(vz_barb/2).^m+1,'','Color',[0.4 0.4 0.4])

plot(vz_bar/2,exp_theodorsen_notiterative,'k')
plot(vz_bar/2,exp_theodorsen_iterative,'g')
plot(vz_bar/2,exp_theodorsen_goldstein,'m')
xlim([0 6])
xlabel('z/D [.]')
ylabel('r/R [.]')
title(sprintf('WakeExpansionDifferentTheoriesCT%2d',CT*100));
legend('Franksen','Rathman','Vortex Cylinder','Theodorsen not iterative','Theodorsen iterative','Theodorsen Goldstein','Location','EastOutside')

