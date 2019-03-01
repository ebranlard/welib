%%
InitClear
require('BEM','v05');
require('WTlib','v06');
require('Wind','v01');

% sWT='SB2';      Format='xblade';
% sWT='Riso10MW'; Format='hawc';
% sWT='NREL5MW_default'; Format='hawc'; 
sWT='NTK500'; Format='hawc'; 

vWS=4:1:25;


[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
WT.Nacelle.tilt=0;
[ WT ]   = fSetRotorGrid(80,WT);
[ Algo ] = fInitAlgo();
Algo.bReInterp=0;
% --------------------------------------------------------------------------------
% --- Computing Power curve based on Turbine specs (vSIMRef), interpolated on vWS
% --------------------------------------------------------------------------------
R=fWTPowerCurve('BEM',WT ,vWS,1,Algo)

%%
figure, fplotCodesComparison('WS','Thrust',{R},'','','',1,1,[],[],'','');
figure, fplotCodesComparison('WS','CT',{R},'','','',1,1,[],[],'','');


%%
