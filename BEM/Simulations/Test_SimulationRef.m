%%
InitClear
setFigurePath('./')
require('BEM','v05');
require('WTlib','v06');
require('Wind','v01');

% sWT='SB2'; Format='xblade'; 
% sWT='Riso10MW'; Format='hawc'; 
sWT='NREL5MW_default'; Format='hawc'; 
% sWT='NTK500'; Format='hawc'; 

[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
WT.Nacelle.tilt=0;
[ WT ]   = fSetRotorGrid(20,WT);
[ Algo ] = fInitAlgo();
% Algo.bRough=0;
Algo.bReInterp=0;


% --------------------------------------------------------------------------------
% --- Finding Turbine specifications 
% --------------------------------------------------------------------------------
Opts.TipSpeedMax=72;
Opts.OmegaMin=6.9*2*pi/60;
Opts.WS_Startup=3;
Opts.vWS_out=3:1:25;
Opts.vLambda=1:1:10;
Opts.vPitch=-2:1:15;
Opts.Pref=WT.Spec.P_rated;
Opts.bOptimBelowPref=true;
Opts.bPlot=1;
[R]= fWTSimRef(WT,Opts,Algo);


% --------------------------------------------------------------------------------
% --- Computing Power curve based on foud specifications
% --------------------------------------------------------------------------------
%%
InitCell
% R2=fWTFindPitch(WT ,4:2:25,Pref,1,Algo);
WT.Spec.vSIMRef=R.vSIMRef;
R2=fWTPowerCurve('BEM',WT ,Opts.vWS_out,1,Algo); %1 for bPlot
figure, fplotCodesComparison('WS','CT',{R2},'','','',1,1,[],[],'','');
%%
% figure,

% figure, fplotCodesComparison('WS','CT',{R2,R},'','','',1,1,[],[],'','');
% 
% %%
% legds={}
%%
% disp('Watch out you are going to clear the workspace')
% pause
% load('WTSimRefWorkspace.mat')
% %%
% figure, fplotCodesComparison('r','CTloc',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),axis ij,box on,colorbar
% figure, fplotCodesComparison('r','a',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),axis ij,box on,colorbar
% figure, fplotCodesComparison('r','Gamma',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),box on,colorbar
%%




