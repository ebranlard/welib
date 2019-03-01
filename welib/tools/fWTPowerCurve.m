function [R]= fWTPowerCurve(Code,WT,vWS,bPlot,Algo)
% 
% Examples:
%  [ WT ]   = fInitWT( 'NTK500' , 'hawc' ,PATH.DATA_WT);
%  [ WT ]   = fSetRotorGrid(30,WT);
%  fWTPowerCurve('BEM',WT ,4:2:24,1,fInitBEMAlgo())

if min(WT.Spec.vSIMRef(:,1))>min(vWS) || max(WT.Spec.vSIMRef(:,1))<max(vWS)
    warning('Not enough input data in SimRef, will have to extrapolate');
end
pitch=interp1(WT.Spec.vSIMRef(:,1), WT.Spec.vSIMRef(:,3),vWS,'cubic','extrap');
rpm=interp1(WT.Spec.vSIMRef(:,1), WT.Spec.vSIMRef(:,2),vWS,'cubic','extrap');
vSIM=[vWS(:) rpm(:) pitch(:)];
[ Sim ]  = fInitSim(WT);
[ Sim ]  = fSetSim( Sim, WT, vSIM ); 
[ Wind ] = fInitWind(  ); 
% keyboard
%% Simulations
[ R ]=fWTSimulation(Code,WT,Sim,Wind,Algo);

if(bPlot)
    Codes={R}; legds={Code};
    colrs=fColrs(1:4);
    sty={'-','+-','--'};
    figure, fplotCodesComparison('WS','PITCH',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','RPM',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('lambda','CP',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('lambda','CT',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','Thrust',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','Edge',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','Flap',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','Power',Codes,legds,colrs,sty,1,1,[],[],'','')
    dispatchFigs(1);
end

