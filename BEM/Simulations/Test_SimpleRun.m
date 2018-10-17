%% Initialization
InitClear;
require('WTlib','v05');
require('BEM','v05');

%% Params
nGrid=30;
% sWT='NTK500p'; Format='hawc'; 
sWT='NTK500';  Format='hawc'; 
% sWT='SB2XLL';  Format='xblade'; 




%% Initialization
[ WT ]   = fInitWT( sWT , Format,PATH.DATA_WT);
[ WT ]   = fSetRotorGrid(nGrid,WT);

lambda=8;
U0=8;
RPM=lambda*U0/WT.Rotor.R*60/2/pi;
[ Sim ]  = fInitSim( WT , [U0 RPM 0]);
[ Wind ] = fInitWind( Sim );

% Setting algorithm parameters
[ Algo ]   = fInitAlgo();
Algo.BEM.bTipLoss=1;

%% Simulation
[ BEM ] = fRunBEM(WT,Sim,Wind,Algo);

%% Plotting results
figure
plot(BEM.r,BEM.Gamma)

figure
plot(BEM.r,BEM.Cl)


dispatchFigs(1)
