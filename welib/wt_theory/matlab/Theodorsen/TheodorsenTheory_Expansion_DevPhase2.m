%% README
% This script is a level Higher than the Expansion_DevPhase0 - It's usefull to investigate iteration, FarwakeParams etc
% For more user friendly script, go to the dedicated WakeExpansion folder and FarWakeParams
% 
% TODO: this script should maybe use fTheodorsenFarWakeParam for the iterative procedure. 
%       It's prbably this script which has been used to derive the iterative method in fTheodorsenFarWakeParam
% The article plot for the torque conference wass done with this script (first part of the script)


%% Init
InitClear
require('OPTIMCIRC','v00');  % v01 uses a mex, use v00 if it doesn't work

setFigurePath({'./figs/' , '/work/publications/articles/2012-tiploss-theoretical/figs/'});
setFigureLatex(1);
R0=100;
U0=10;
nB=3;
nr=100;
ntheta=3000;
z_inf=100;
bWT=1;

%% Standalone for Ewan test
addpath('/work/lib/WTTheory/Expansion/'); 

lambda=7.06;
CT=0.737;
nB=3;
R0=27;
U0=10
a=1/2*(1-sqrt(1-CT));
 [ Expansion vz_bar R_over_RW_Wald] = fTheodorsenExpansionWrap(lambda,nB,U0,R0,a,CT );
 [ Expansion2 vz_bar2] = fExpansionVortexRings(CT );

figure,hold all
plot(vz_bar/2,Expansion)
plot(vz_bar2/2,Expansion2)
xlim([0 5])

%% DONE AS A FUNCTION OF just l_bar, ct or cs not taken into account -> OF COURSE NOT ITERATIVE THEN
vlbar=1./[10 4 2 1];
w_bar=0.4;

legds={};
for il=1:length(vlbar)
    l_bar=vlbar(il);
    vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
    [ K kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );
    [ vz_bar ExpFactor Y Y1 Y2 kappa] = fTheodorsenExpansion(l_bar,vx,K,nB,z_inf,ntheta,bWT );

    figure(3),set(3,'DefaultAxesColorOrder',repmat(linspace(0,0.7,length(vlbar)),3,1)');
%     hold all, grid on,box on,   plot(vz_bar,(1+ExpFactor)/(1+max(ExpFactor))),title('Expansion');
    hold all, grid on,box on,   plot(vz_bar/2,1+cs*ExpFactor),title('Expansion');
    xlim([0 10])
    ylim([1 1.3])
    [l_bar w_bar cs]

    legds{end+1}=sprintf('$1/\\bar{l}=%d$ - $c_t=%.2f$',1/vlbar(il),cs);
end
figure(3)
legend(legds)
title('WakeExpansionTheodorsenLBAR')
xlabel('$z/D$ [.]')
ylabel('$R_w/R$ [.]')
dispatchFigs(1)


%% DONE AS A FUNCTION OF FAR WAKE PARAMS -> OF COURSE NOT ITERATIVE THEN
vlbar=1./[8 8 3 3];
vwbar=[0.6 0.4 0.4 0.2];

legds={};
for il=1:length(vlbar)
    l_bar=vlbar(il);
    w_bar=vwbar(il);

    vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
    [ K kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );
    [ vz_bar ExpFactor Y Y1 Y2 kappa] = fTheodorsenExpansion(l_bar,vx,K,nB,z_inf,ntheta,bWT );

    figure(3),set(3,'DefaultAxesColorOrder',repmat(linspace(0,0.7,length(vlbar)),3,1)');
    hold all, grid on,box on,   plot(vz_bar/2,1+cs*ExpFactor),title('Expansion');
    xlim([0 10])
    ylim([1 1.5])
    [l_bar w_bar cs]
    legds{end+1}=sprintf('$1/\\bar{l}=%d$ - $\\bar{w}=%0.2f$ - $c_t=%.2f$',1/vlbar(il),vwbar(il),cs);
end


%
figure(3)
legend(legds)
title('WakeExpansionTheodorsenFarWakeParams')
xlabel('$z/D$ [.]')
ylabel('$R_w/R$ [.]')

%% HERE THIS IS DONE WITH A WEIRD RELATION WITH RESPECT TO LAMBDA
vlambda=[12 9 7 5 3];
vCT=[1.1 0.8 0.6 0.3 0.1];
va=[0.5 0.33 0.25 0.12 0.05];

legds={};
for il=1:length(vlambda)
    U=U0*(1-va(il));
    l_bar=(1-va(il))/vlambda(il);
    CT=vCT(il);
    w_bar=2*va(il);

    vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 
    [ K kappa eps_over_k epsz cs cp e eta a0 ] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );

    [ vz_bar ExpFactor Y Y1 Y2 kappa] = fTheodorsenExpansion(l_bar,vx,K,nB,z_inf,ntheta,bWT );

    Rw=ExpFactor(end)*vCT(il)*R0+R0;

    figure(1),hold all, grid on,box on,   plot(vx,K),title('K')
    figure(2),hold all, grid on,box on,   plot(vz_bar,ExpFactor),title('ExpansionFactor');
    xlim([0 10])
    figure(3),set(3,'DefaultAxesColorOrder',repmat(linspace(0,0.7,length(vlambda)),3,1)');
    hold all, grid on,box on,   plot(vz_bar/2,1+ExpFactor*vCT(il)),title('Expansion');
    xlim([0 10])

    [l_bar w_bar CT cs Rw]

    legds{end+1}=sprintf('\\lambda=%d - C_T=%0.2f',vlambda(il),vCT(il));
    dispatchFigs(1)
end

dispatchFigs(1)

%
close(1:2)
figure(3)
legend(legds)
title('WakeExpansionTheodorsenLambdaCTNotIterative')




%%
[ Expansion vz_bar] = fTheodorsenExpansionWrap(vlambda(1),nB,U0,R0,va(1),vCT(1));
figure(4),hold all, grid on,   plot(vz_bar/2,Expansion),title('Expansion');



%% The smae but iterative !!!!!!!!!!!!!!!!!!!!!!!!!!!! to be reviewed
% This has been used to establish fTheodorsenFarWakeParams

vlambda=[12 9 7 5 3];
vCT=[1.1 0.8 0.6 0.3 0.1];
va=[0.5 0.33 0.25 0.12 0.05];

vlambda=[12];
vCT=[1.1];
va=[0.5];


vlambda=[12];
vCT=[1.1];
va=[0.5];



legds={};
for il=1:length(vlambda)
    U=U0*(1-va(il));
    l_bar=(1-va(il))/vlambda(il);
    lambda=vlambda(il);
    CT=vCT(il);
    w_bar=2*va(il);
    
    Rw=1.1*R0; % that's the approximation for first round
   
    cs=CT*R0^2/Rw^2;

    while abs(cs-CT)/max(cs,CT)>10^-6
        vx=linspace(0,1,nr);  % Note: there is a singularity at x=1, see the denominator of y1 for theta=tau=0. The denaminator cancels out. Normal, vortex singularity 

        
        l_bar=(1-w_bar/2)/lambda;  % that's probably the link given by Okulov and Sorensen.

        [ K kappa eps_over_k epsz cs] = fTheodorsenAtOperatingPoint( l_bar, nB, vx ,w_bar );
        [ vz_bar ExpFactor ] = fTheodorsenExpansion(l_bar,vx,K,nB,z_inf,ntheta,bWT );


        Rw=ExpFactor(end)*cs*R0+R0;
        % for next round
        cs_new=CT*R0^2/Rw^2
        l_bar_new=l_bar
        l_bar_sor=(1-w_bar/2)/lambda;  % that's probably the link given by Okulov and Sorensen.

        figure(2),hold all, grid on,   plot(vz_bar/2,ExpFactor),title('ExpansionFactor');
        xlim([0 10])
        figure(3),hold all, grid on,   plot(vz_bar/2,1+ExpFactor*cs),title('Expansion');
        xlim([0 10])
        legds{end+1}=sprintf('\\lambda=%d - C_T=%0.2f',vlambda(il),vCT(il));
        dispatchFigs(1)

        [l_bar w_bar CT cs Rw]
    end
end
legend(legds)
dispatchFigs(1)
%%
