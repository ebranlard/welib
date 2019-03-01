%% Initialization
InitClear;
require('WTlib','v05');
require('BEM','v05');
require('OPTIMCIRC','v01');  
PATH.CirculationFamily='/mnt/DataWin/Work/2011-Siemens/code/CirculationFamilyCurves/';
require('CirculationFamily','');  
sWT='SB1'; sFormat='xblade';
% WT init
R=1;  % probably does not work if R !=1
[ WT ]   = fInitWT( sWT ,sFormat,PATH.DATA_WT);
WT.Rotor.cone=0;
WT.Nacelle.tilt=0;
WT.Sources.Rotor.chord=WT.Sources.Rotor.chord/WT.Rotor.R;
WT.Sources.Rotor.twist=WT.Sources.Rotor.twist*0+20;
WT.Rotor.R=R; WT.Rotor.rhub=0.00; WT.Rotor.BladeLength=R; WT.Rotor.R_coned=R; WT.Rotor.SweptArea=pi*R^2;
WT.Sources.Rotor.r=linspace(WT.Rotor.rhub,R,length(WT.Sources.Rotor.r));
Scale=60/WT.Rotor.R;
[ WT ] = fSetRotorGrid(20,WT);
% BEM Algorithm init
[ Algo ] = fInitAlgo();
Algo.bPrescribedGamma=1;
% Simulation/Circulation Params
x0=0.1; x2=0; y3=0; t0=0.3;
lambda=3;
x0=0.2; x2=0.5; y3=0.5; t0=0.2;
lambda=3;
U0=8;
RPM=lambda*U0/WT.Rotor.R*60/2/pi;
[ Sim ]  = fInitSim( WT , [U0 RPM -10]);
[ Wind ] = fInitWind( Sim );


%% Prescribing circulation - Amplitude
rGam0=linspace(0,R,100);
% My parametrization
Amplitude=2;
[sse Gam0]=fFitGamma([x0 x2 y3 t0],rGam0,rGam0*0);
Algo.r_PrescribedGamma=rGam0;
Algo.PrescribedGamma=Gam0*Amplitude;
%Goldstein
CT=0.5;
nB=WT.Rotor.nB;
[w_bar l_bar CT G]=fGoldsteinFarWakeParams(CT,lambda,nB);
Gamma=G/nB*w_bar*U0*2*pi*R*l_bar;
Algo.r_PrescribedGamma=rGam0;
Algo.PrescribedGamma=Gamma;


figure,hold all
plot(Algo.r_PrescribedGamma,Algo.PrescribedGamma);
Sim.Name=sprintf('TipLossDB_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_',x0,x2,y3,t0,lambda,Amplitude);
 Algo.BEM.bTipLoss=0;
[ BEM ] = fRunBEM(WT,Sim,Wind,Algo);
BEM.CT
plot(BEM.r,BEM.Gamma);
plot(BEM.r,BEM.a);
plot(BEM.r,BEM.aprime);
plot(BEM.r,BEM.CTloc,'k');


%%
