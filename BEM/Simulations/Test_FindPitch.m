%%
InitClear
require('BEM','v05');
require('WTlib','v05');

sWT='SB2'; Format='xblade'; 
sWT='Riso10MW'; Format='hawc'; 
% sWT='NTK500'; Format='hawc';
[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
WT.Nacelle.tilt=0;
[ WT ]   = fSetRotorGrid(12,WT);
[ Algo ] = fInitAlgo();
Algo.bReInterp=0;

R=fWTFindPitch(WT ,4:1:25,WT.Spec.P_rated,1,Algo);
WT.Spec.vSIMRef=R.SIMRef;
fWTPowerCurve('BEM',WT ,4:1:25,1,Algo)

