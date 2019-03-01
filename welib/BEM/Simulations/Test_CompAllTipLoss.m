%% Initialization - Typical constants
InitClear;
require('WTlib','v05');
require('BEM','v05');
require('OPTIMCIRC','v01');
setFigureLatex(1);
global legds;
%%
sWT='NRELShen'; Format='hawc';
sWT='SB2';  Format= 'xblade';

[ WT ]   = fInitWT( sWT, Format, PATH.DATA_WT);
 WT.Rotor.cone=0;% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% case definition
if isequal(sWT,'NRELShen')
    i=0; load('data/Shen6.mat'); load('data/NrelExp6.mat'); Nrel=NrelExp6; Shen=Shen6;
% i=1; load('data/Shen10.mat'); load('data/NrelExp10.mat'); Nrel=NrelExp10;  Shen=Shen10;
vSIM=WT.Spec.vSIM(WT.Spec.iSimDef+i,:);
vSIM(3)=4;
else
vSIM=WT.Spec.vSIM(WT.Spec.iSimDef,:);
end
[ Sim ]  = fInitSim( WT , vSIM );
[ Wind ] = fInitWind( Sim );
[ Algo ] = fInitAlgo();
[ WT ] = fSetRotorGrid(20,WT);



Rotor=WT.Rotor;
Algo.TIDrag=1;
Algo.BEM.CTCorrection='GlauertCT';
if isequal(sWT,'SB2') 
    Algo.bRoughProfiles=1;
else
    Algo.bReInterp=0;
    Algo.RoughProfiles=0;
end
XLIM=[0 0.99];

%%
disp('No loss')
Algo.BEM.bTipLoss=0;
[ BEMNoLoss ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('Glauert')
Algo.BEM.bTipLoss=1;
Algo.BEM.TipLossMethod='Glauert';
[ BEMGl ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('Xu Sankar')
Algo.BEM.TipLossMethod='XuSankar';
[ BEMXS ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('Prandtl')
Algo.BEM.TipLossMethod='Prandtl';
[ BEMPr ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('Lindenburg')
Algo.BEM.TipLossMethod='Lindenburg';
[ BEMLi ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('GoldsteinSimple')
Algo.BEM.TipLossMethod='GoldsteinSimple';
[ BEMGo ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('Shen')
Algo.bTipLossCl=1;
Algo.BEM.TipLossMethod='Shen';
[ BEMSh ] = fRunBEM(WT,Sim,Wind,Algo);

%%
disp('TipLossDB')
Algo.bTipLossCl=0;
Algo.BEM.TipLossMethod='TipLossDB';
[ BEM BEMManu ] = fRunBEM(WT,Sim,Wind,Algo);








% --------------------------------------------------------------------------------
% --- Plots using bPlotCodeCompBEM
% --------------------------------------------------------------------------------
%% Tip-Loss
figure,hold all,grid,box,legds={};
fPlotCodeCompBEM('F',BEMGl,BEMNoLoss,BEMManu,BEMSh,BEMXS,WT.Rotor.R)

legend(legds,3)
ylabel('Tip-Loss factor $F$ [.]')
xlabel('$r/R$ [.]')
xlim(XLIM)
title(sprintf('%s_F',Sim.Name))

%%
figure,hold all,grid,box,legds={};
fPlotCodeCompBEM('Re',BEMGl,BEMNoLoss,BEMManu,BEMSh,BEMXS,WT.Rotor.R)
legend(legds,3)
ylabel('Re [.]')
xlabel('$r/R$ [.]')
xlim(XLIM)
title(sprintf('%s_Re',Sim.Name))

%%
figure,hold all,grid,box,legds={};
fPlotCodeCompBEM('alpha',BEMGl,BEMNoLoss,BEMManu,BEMSh,BEMXS,WT.Rotor.R)
legend(legds,3)
ylabel('alpha [.]')
xlabel('$r/R$ [.]')
xlim(XLIM)
title(sprintf('%s_alpha',Sim.Name))
%%
figure,hold all,grid,box,legds={};
fPlotCodeCompBEM('Gamma',BEMGl,BEMNoLoss,BEMManu,BEMSh,BEMXS,WT.Rotor.R)
legend(legds,3)
ylabel('Gamma [.]')
xlabel('$r/R$ [.]')
xlim(XLIM)
title(sprintf('%s_alpha',Sim.Name))


%%
figure,hold all,grid,box,legds={};
fPlotCodeCompBEM('Cl',BEMGl,BEMNoLoss,BEMManu,BEMSh,BEMXS,WT.Rotor.R)
legend(legds,3)
ylabel('Cl [.]')
xlabel('$r/R$ [.]')
xlim(XLIM)
title(sprintf('%s_Cl',Sim.Name))
%%
figure,hold all,grid,box,legds={};
fPlotCodeCompBEM('Pn',BEMGl,BEMNoLoss,BEMManu,BEMSh,BEMXS,WT.Rotor.R)
ylabel('$P_n$ [N/m]')
xlabel('$r/R$ [.]')
xlim(XLIM)
title(sprintf('%s_alpha',Sim.Name))
legds{end+1}='Exp';
legend(legds,3)
if isequal(sWT,'NRELShen')
    plot(Nrel(:,1),-Nrel(:,2),'s')
end
grid on
box on




% --------------------------------------------------------------------------------
% --- Manual Plots
% --------------------------------------------------------------------------------
%set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1], 'DefaultAxesLineStyleOrder','-|--|:|-.')
figure
hold all
plot(Rotor.r/Rotor.R,-BEMGl.Pn,'k')
plot(Rotor.r/Rotor.R,-BEMNoLoss.Pn,'k--')
plot(Rotor.r/Rotor.R,-BEMGo.Pn,'-','Color',fColrs(1))
plot(Rotor.r/Rotor.R,-BEMPr.Pn,'--','Color',fColrs(1))
plot(Rotor.r/Rotor.R,-BEMLi.Pn,'-','Color',fColrs(2))
plot(Rotor.r/Rotor.R,-BEMSh.Pn,'-','Color',fColrs(3))
plot(Rotor.r/Rotor.R,-BEMXS.Pn,'--','Color',fColrs(3))
plot(Rotor.r/Rotor.R,-BEMManu.Pn,'-','Color',fColrs(1),'LineWidth',2)
xlim([0.7 1])
grid on
box on
xlabel('$r/R$ [.]')
ylabel('Normal load $P_n$ [.]')
title('AllTipLossesPn')
legend('Glauert','No Loss','Goldstein','Prandtl','Lindenburg','Shen','Xu Sankar',0)


%%
%set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1], 'DefaultAxesLineStyleOrder','-|--|:|-.')
figure
hold all
%plot(Rotor.r/Rotor.R,BEMNoLoss.F,'k--')
plot(Rotor.r/Rotor.R,BEMGl.F,'k')
plot(Rotor.r/Rotor.R,BEMGo.F,'-','Color',fColrs(1))
plot(Rotor.r/Rotor.R,BEMPr.F,'--','Color',fColrs(1))
% plot(Rotor.r/Rotor.R,BEMSh.F,':','Color',fColrs(2))
plot(Rotor.r/Rotor.R,BEMLi.F,'-','Color',fColrs(2))
plot(Rotor.r/Rotor.R,BEMXS.F,'--','Color',fColrs(3))
plot(Rotor.r/Rotor.R,BEMManu.F,'-','Color',fColrs(1),'LineWidth',2)
title('AllTipLossesF')
grid on
box on
xlabel('$r/R$ [.]')
ylabel('Tip-loss factor $F$ [.]')
legend('Glauert','Goldstein','Prandtl','Lindenburg','Xu Sankar',0)
%%

figure
hold all
plot(Rotor.r/Rotor.R,BEMGl.a,'k')
plot(Rotor.r/Rotor.R,BEMNoLoss.a,'k--')
plot(Rotor.r/Rotor.R,BEMGo.a,'-','Color',fColrs(1))
plot(Rotor.r/Rotor.R,BEMPr.a,'--','Color',fColrs(1))
plot(Rotor.r/Rotor.R,BEMLi.a,'-','Color',fColrs(2))
plot(Rotor.r/Rotor.R,BEMSh.a,'-','Color',fColrs(3))
plot(Rotor.r/Rotor.R,BEMXS.a,'--','Color',fColrs(3))
plot(Rotor.r/Rotor.R,BEMManu.a,'-','Color',fColrs(1),'LineWidth',2)
ylim([0 1])
xlim([0.7 1])
title('AllTipLossesA')
grid on
box on
xlabel('$r/R$ [.]')
ylabel('Axial induction $a_B$ [.]')
legend('Glauert','No Loss','Goldstein','Prandtl','Lindenburg','Shen','Xu Sankar',0)


%%

figure
hold all
plot(Rotor.r/Rotor.R,BEMGl.alpha,'k')
plot(Rotor.r/Rotor.R,BEMNoLoss.alpha,'k--')
plot(Rotor.r/Rotor.R,BEMGo.alpha,'-','Color',fColrs(1))
plot(Rotor.r/Rotor.R,BEMPr.alpha,'--','Color',fColrs(1))
plot(Rotor.r/Rotor.R,BEMLi.alpha,'-','Color',fColrs(2))
plot(Rotor.r/Rotor.R,BEMSh.alpha,'-','Color',fColrs(3))
plot(Rotor.r/Rotor.R,BEMXS.alpha,'--','Color',fColrs(3))
% ylim([0 1])
xlim([0.7 1])
title('AllTipLossesAlpha')
grid on
box on
xlabel('$r/R$ [.]')
ylabel('Angle of attack $\alpha$ [.]')
legend('Glauert','No Loss','Goldstein','Prandtl','Lindenburg','Shen','Xu Sankar',0)

%%
CPS=[BEMNoLoss.CP  BEMGl.CP  BEMGo.CP  BEMPr.CP  BEMLi.CP BEMSh.CP  BEMXS.CP];
mean=sum(CPS)/7;
M=max(abs(mean-CPS));
reldiff=(CPS-mean)/M;

%%
figure
hold all
for i=1:length(CPS)
    bar(i,CPS(i),'FaceColor',fColorGrad(reldiff(i)))
end

% bar([BEMNoLoss.CP BEMPr.CP BEMGo.CP BEMGl.CP BEMSh.CP BEMLi.CP BEMXS.CP],'FaceColor',[fColrs(1);fColrs(2)])
% title('AllTipLossesCp')
set(gca,'XTickLabel',{'','No Loss','Glau.','Gold.','Pran.','Lind.','Shen','Xu\&Sa.',''});
box on
ylabel('Power Coefficient $C_P$ [.]')
title('AllTipLossesCP')
%%

figure
hold all
for i=1:length(CPS)
    bar(i,(CPS(i)-mean)/mean*100,'FaceColor',fColorGrad(reldiff(i)))
end
% bar([BEMNoLoss.CP BEMPr.CP BEMGo.CP BEMGl.CP BEMSh.CP BEMLi.CP BEMXS.CP],'FaceColor',[fColrs(1);fColrs(2)])
% title('AllTipLossesCp')
set(gca,'XTickLabel',{'','No Loss','Glau.','Gold.','Pran.','Lind.','Shen','Xu\&Sa.',''});
box on
ylabel('Relative difference in $C_P$ [\%]')
title('AllTipLossesCPRel')



dispatchFigs(1)
