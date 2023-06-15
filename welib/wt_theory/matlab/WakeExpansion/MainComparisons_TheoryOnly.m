% This script looks only at the wake expansion (it doens' look at tip-losses or BEM etc.)
% There are some nasty copy paste to look at different CT

%% 
InitClear;
PATH.EXPANSION='/work/lib/WTTheory/WakeExpansion/';
require('EXPANSION','');

require('THEODORSEN','v01');
require('OPTIMCIRC','v01_nomex');


setFigurePath({'./' , '/work/publications/articles/2012-tiploss-theoretical/figs/'});

vCT=[0.35 0.65 0.8];
vCT=0.737;
for CT=vCT
    vz_bar=linspace(0,30,100); %z/R !!! not z/D

    Algo.nr=100;
    Algo.ntheta=length(vz_bar); % dangerous
    Algo.z_inf=max(vz_bar);     % dangerous
    Algo.bWT=1;
%     lambda=15;
    lambda=7.06;
    nB=3;
    U0=10;
    R=1;

    [ exp_franksen ] = fExpansionFranksen(CT,0.7,2,vz_bar);
    [ exp_rathmann ] = fExpansionRathmann(CT,vz_bar);
    [ exp_theodorsen ] = fTheodorsenExpansionWrap(lambda,nB,U0,R,1/2*(1-sqrt(1-CT)),CT, Algo); % !!!! mean a, that's not really good
    [ exp_vortexrings ] = fExpansionVortexRings(CT,vz_bar); 


    %%
    figure,hold all,grid on,box on
    plot(vz_bar/2,exp_franksen,'')
    plot(vz_bar/2,exp_rathmann,'')
    plot(vz_bar/2,exp_theodorsen,'k')
    plot(vz_bar/2,exp_vortexrings,'k--')
    k=0.05;m=1;
    plot(vz_bar/2,k*(vz_bar/2).^m+1,'','Color',[0.4 0.4 0.4])
    xlim([0 6])
    xlabel('z/D [.]')
    ylabel('r/R [.]')
    title(sprintf('WakeExpansionDifferentTheoriesCT%2d',CT*100));
    legend('Franksen','Rathman','Theodorsen','Vortex Cylinder','k x^m')

end
%%
