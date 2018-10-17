%% Initialization for the BEM code
InitClear;
require('WTlib','v05');
require('BEM','v05');


%% Params
nGrid=30;
vWS=4:1:24;
% sWT='NTK500p'; Format='hawc';  load('/work/lib/WTlib/v02/MainFindPitch/NTK500p_stall.mat')
sWT='NTK500';  Format='hawc';  load('/work/lib/WTlib/v02/MainFindPitch/NTK500_stall.mat')



%% Initialization
[ WT ]   = fInitWT( sWT , Format,PATH.DATA_WT);
WT.Nacelle.tilt=0;
[ WT ]   = fSetRotorGrid(nGrid,WT);
[ Sim ]  = fInitSim( WT ); %load default specs as sim
[ Sim ]  = fSetSim( Sim, WT, vWS ,WT.Spec.Omega_rated*60/(2*pi), 0,0  ); 
[ Wind ] = fInitWind(  ); 

% Setting algorithm parameters
[ Algo ]   = fInitAlgo();
Algo.bReInterp=0;
Algo.aTol=10^-6;
Algo.relaxation=0.2;
Algo.BEM.bHubLoss=0;
Algo.BEM.CTCorrection='Hawc';



%% Simulation
[ BEM ]=fWTSimulation('BEM',WT,Sim,Wind,Algo);


%% Plotting results
Codes={BEM}; legds={'BEM'};
colrs=fColrs(1:4);
sty={'-','+-','--'};
R=WT.Rotor.R;
figure, fplotCodesComparison('WS','Power',Codes,legds,colrs,sty,1,1,[],[],'','')
plot(par.Uvec,par.P/1000,'k-+');

figure, fplotCodesComparison('WS','CT',Codes,legds,colrs,sty,1,1,[],[],'','')
figure, fplotCodesComparison('lambda','CP',Codes,legds,colrs,sty,1,1,[],[],'','')
figure, fplotCodesComparison('lambda','CT',Codes,legds,colrs,sty,1,1,[],[],'','')
figure, fplotCodesComparison('WS','Thrust',Codes,legds,colrs,sty,1,1,[],[],'','')
figure, fplotCodesComparison('WS','Edge',Codes,legds,colrs,sty,1,1,[],[],'','')
figure, fplotCodesComparison('WS','Flap',Codes,legds,colrs,sty,1,1,[],[],'','')


CodesAll=squeeze(BEM.Results);
figure, fplotCodesComparison('r','Cn',CodesAll,legds,colrs,sty,'scaleByMax',1,[],[],'','',vWS,'U [m/s]')


dispatchFigs(1)

