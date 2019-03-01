%% Initialization - Typical constants
InitClear;
% addpath([MAINPATH 'code/WTlib/']);
% addpath([MAINPATH 'code/BEM/']);
% addpath(CCODEPATH);
require('BEM','v02')
require('WTlib','v02')
require('OPTIMCIRC','v-1')
setFigureLatex(1)
%%
%sGeometry='flatplate';
sGeometry='elliptic';
sGeometry='Manu';
% sGeometry='FlatTUD'; % works better with no viscous model
sGeometry='NB4';
sGeometry='NRELShen';
%  sGeometry='B49';

[ WT ]   = fInitWT( sGeometry ,'xblade',PATH.DATA_WT);
WT.Rotor.cone=0;% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% case definition

[ Algo ] = fInitBEMAlgo();
[ WT ] = fSetRotorGrid(60,WT);
Rotor=WT.Rotor;


PITCH=WT.Spec.vSIM(WT.Spec.iSimDef,3);
RPM=WT.Spec.vSIM(WT.Spec.iSimDef,2);
vWS=3:0.5:25;


for i=1:length(vWS)
    [ Sim ]  = fInitSim( WT , [vWS(i) RPM PITCH] );
    [ Wind ] = fInitWind( Sim );
    %
    Algo.TipLossMethod='Glauert';
    Algo.bTipLoss=0;
    [ BEMNoLoss(i) ] = fRunBEM(WT,Sim,Wind,Algo);

    %
    Algo.bTipLoss=1;


    Algo.TipLossMethod='GoldsteinSimple';
    [ BEMGo(i) ] = fRunBEM(WT,Sim,Wind,Algo);

    %
    Algo.TipLossMethod='Manu';
    [ ~, BEMMa(i) ] = fRunBEM(WT,Sim,Wind,Algo);
    %
    Algo.TipLossMethod='Glauert';
    [ BEMGl(i) ] = fRunBEM(WT,Sim,Wind,Algo);


    Algo.TipLossMethod='Shen';
    [ BEMSh(i) ] = fRunBEM(WT,Sim,Wind,Algo);

    %
    Algo.TipLossMethod='XuSankar';
    [ BEMXS(i) ] = fRunBEM(WT,Sim,Wind,Algo);

    Algo.TipLossMethod='Prandtl';
    [ BEMPr(i) ] = fRunBEM(WT,Sim,Wind,Algo);
    Algo.TipLossMethod='Lindenburg';
    [ BEMLi(i) ] = fRunBEM(WT,Sim,Wind,Algo);


end



save('data/AllTipLossComparisonAEP_Nov2012.mat')
%%

%% Some definitions
% load('data/AllTipLossComparisonAEP.mat')
Power{1}=[BEMGl(:).Power];
Power{2}=[BEMNoLoss(:).Power];
Power{3}=[BEMGo(:).Power];
Power{4}=[BEMPr(:).Power];
Power{5}=[BEMLi(:).Power];
Power{6}=[BEMSh(:).Power];
Power{7}=[BEMXS(:).Power];
Power{8}=[BEMMa(:).Power];
WS=vWS;
%Weibull param
%figure 1
A2=6.8;  %6.8
k2=2.16;

%figure 2
A1=8.9;  %6.8
k1=2.16;


p = mywblpdf(WS,A1,k1);
nHpY=24*365.25;% Hours per year
vSty={'-','--','-','--','-','-','--','-'};
vColrs(1,:)=[0 0 0];
vColrs(2,:)=[0 0 0];
vColrs(3,:)=fColrs(1);
vColrs(4,:)=fColrs(1);
vColrs(5,:)=fColrs(2);
vColrs(6,:)=fColrs(3);
vColrs(7,:)=fColrs(3);
vColrs(8,:)=fColrs(4);



%% Weibull
WS2=linspace(0,25,50);

p1 = mywblpdf(WS2,A1,k1);
p2 = mywblpdf(WS2,A2,k2);
figure,clf,hold all,box on,grid on
%bar(WS,p*nHpY,'FaceColor',fColrs(3))
plot(WS2,p1*nHpY,'Color',fColrs(1))
plot(WS2,p2*nHpY,'Color',fColrs(2))
xlim([0 25])
xlabel('Wind Speed [m/s]')
ylabel('Hours of operation per year [h]')
title('AEPWindContent')
legend('A=8.9 - k=2.16','A=6.8 - k=2.16',0)




%%
% load('data/NRELExp.mat')

figure
hold all
for i=1:length(Power)
    plot(vWS,Power{i}/1000,vSty{i},'Color',vColrs(i,:))
end
ylim([0 13])
grid on
box on
xlabel('$r/R$ [.]')
ylabel('Aerodynami Power $P$ [kW]')
title('AllTipLossesPower')
legend('Glauert','No Loss','Goldstein','Prandtl','Lindenburg','Shen','Xu Sankar','Manu','Location','SouthEast')



%%


% figure,clf,hold all,box on,grid on
% %bar(WS,p*nHpY,'FaceColor',fColrs(3))
% plot(WS,p*nHpY)
% xlim([0 25])
% xlabel('Wind Speed [m/s]')
% ylabel('Hours of operation per year [h]')
% title('AEPWindContent')

%%
% figure(2),clf,hold all,box on,grid on
eta_elec=0.97;
eta_oper=0.99;
for i=1:length(Power)
    P=Power{i}*eta_elec;
    WS2=0:0.2:26;
    Pinterp=interp1(WS,P,WS2,'spline','extrap');
    I=find(WS2<4 | WS2>25);
    Pinterp(I)=0;
    p1 = mywblpdf(WS2,A1,k1);
    p2 = mywblpdf(WS2,A2,k2);
    MWh1(i,:)=[p1.*Pinterp*nHpY/10^6]*eta_oper;
    MWh2(i,:)=[p2.*Pinterp*nHpY/10^6]*eta_oper;
    AEP1(i)=trapz(WS2,p1.*Pinterp*nHpY)/10^6*eta_oper;
    AEP2(i)=trapz(WS2,p2.*Pinterp*nHpY)/10^6*eta_oper;
    %     plot(vWS,MWh(i,:)/1000,vSty{i},'Color',vColrs(i,:))
end

%%
flag=1;
%


if flag==1
    AEP=AEP1;
else
    AEP=AEP2;
end



meanAEP=mean(AEP);
refVal=meanAEP;  % what I used to do
% refVal=AEP(1);  % no loss values

M=max(abs(refVal-AEP));
reldiff=(AEP-refVal)/M;

figure, hold all, box on
I=[2 1 3:length(AEP)];
for i=1:length(I)
    bar(i,AEP(I(i)),'FaceColor',fColorGrad(reldiff(I(i))))
end
set(gca,'XTickLabel',{'','No Loss','Glau.','Gold.','Pran.','Lind.','Shen','Xu\&Sa.','Bran.',''});
ylabel('AEP [MWh/year]')
title(['AllTipLossesAEPnew' num2str(flag)])

figure, hold all, box on
for i=1:length(I)
    bar(i,(AEP(I(i))-refVal)/refVal*100,'FaceColor',fColorGrad(reldiff(I(i))))
end
set(gca,'XTickLabel',{'','No Loss','Glau.','Gold.','Pran.','Lind.','Shen','Xu\&Sa.','Bran.',''});
ylabel('Relative difference in AEP [\%]')
title(['AllTipLossesAEPRelnew' num2str(flag)])
