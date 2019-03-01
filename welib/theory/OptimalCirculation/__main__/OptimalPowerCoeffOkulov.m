InitClear;
setFigurePath({'./' , '/work/publications/articles/2012-tiploss-theoretical/figs/'});
setFigureTitle(0);
setMatFigure(1);
setFigureLatex(1);
% 
% Fine this script loops on lambda and uses the link between lambda lbar and w_bar.
% It would have been better just to do it as function of l_bar, as the first article from 2007 "optimal operational regime".
% Also I could have used fzero or fminsearch for the max cp to find w_bar. This 2/3 puzzles me. 

%% FLAG 
bPlot=1; % high resolution for plot, low res for table
i3=2; % three blades, at index 2

VB=[1 3 5 100];
VN=[100 100 100 100 100];



%% TEST new function
lbar=1/2;
nB=12;
N=100;
Vr=linspace(0,1,N);
tic()
Gref=fGoldsteinFactor( lbar, nB, Vr );
toc()
tic()
G=fGoldsteinFactor_Matlab( lbar, nB, Vr );
toc()
% tic()
% G2=fGoldsteinFactor_MatlabOld( lbar, nB, Vr );
% toc()

close all
figure,hold all
plot(Vr,Gref,'b-')
plot(Vr,G,'k+')
% plot(Vr,G2,'go')

%% PART 1 : As function of l_bar - No relation with near-wake parameters - First article Okulov 2007 Optimal
tic();
%%
close all

Vl_bar_inv=0.5:0.5:15;
vl_bar=1./Vl_bar_inv/(2*pi);
%VN=[100 100 70 30];
OkulovVw=zeros(length(VB),length(vl_bar));
OkulovCp=zeros(length(VB),length(vl_bar));
OkulovCT=zeros(length(VB),length(vl_bar));
s=-1; % s=-1: WT - s=+1: Propellers
for iB=1:length(VB)
    B=VB(iB);
    Vr=linspace(0,1,VN(iB)); 

    for il=1:length(Vl_bar_inv)
        l_bar=vl_bar(il);        
        G=fGoldsteinFactor( l_bar,B,Vr );
        I1=2*trapz(Vr,G.*Vr);
        I2=2*trapz(Vr,G*l_bar^2.*Vr./(l_bar^2+Vr.^2));  % Corrected version of Okulov's equation
%         I2=2*trapz(Vr,G*l_bar.*Vr.^2./(l_bar^2+Vr.^2)); OKULOV VERSION 
        I3=2*trapz(Vr,G.*Vr.^3./(l_bar^2+Vr.^2));
        if(s==1)
            w_bar=1/(3*I3)*(-I1-I3+sqrt(I1^2+I3^2-I1*I3));
        else
            w_bar=1/(3*I3)*(I1+I3-sqrt(I1^2+I3^2-I1*I3));
        end
        a=w_bar*mean(G); % should be replaced by trapz
        sigma=(1+s*a)/(1+2*s*a); % Correction version of Okulov's
%         w_bar=1/(3*I3)*(sqrt(I2^2+I1*I3)-I1-I3); % OKULOV VERSION
        fprintf('B=%03d - lbar=%04.1f - wbar %.3f\n',B,l_bar,w_bar);
        OkulovVw(iB,il)=w_bar;
        OkulovVa(iB,il)=a;
        OkulovVsigma(iB,il)=sigma;
%       Okulov  sigma=1;
        OkulovCp(iB,il)=2*sigma*w_bar*(1+s*w_bar)*(I1+s*w_bar*I3);
%       Okulov  Cp(iB,il)=2*sigma*w_bar*(1+s*w_bar)*(I1+s*w_bar*I3); %OKULOV VERSION
        OkulovCT(iB,il)=2*sigma*(w_bar*(1+s*w_bar)*I1-s*w_bar^2*I2);
    end
end
legds={};
for iB=1:length(VB)
    legds{end+1}=sprintf('B=%d',VB(iB));
    figure(1)
    hold all
    plot([0 Vl_bar_inv],[0 OkulovCp(iB,:)])
    
    figure(2)
    hold all
    plot([0 Vl_bar_inv],[0 OkulovCT(iB,:)])
%     
    figure(3)
    hold all
    plot([0 Vl_bar_inv],[0 OkulovVa(iB,:)])
% 
%     
    figure(4)
    hold all
    plot([0 Vl_bar_inv],[1 OkulovVsigma(iB,:)])
end


%% PART 2: Okulov-Sorensen lambda - It has to be iterative since l_bar depends on w_bar
if bPlot
    Vlambda=0.2:0.2:15;
else
    Vlambda=1:1:15;
end
Vw=zeros(length(VB),length(Vlambda));
Cp=zeros(length(VB),length(Vlambda));
CT=zeros(length(VB),length(Vlambda));
for iB=1:length(VB)
    B=VB(iB);
    Vr=linspace(0,1,VN(iB));  
    for il=1:length(Vlambda)
        lambda=Vlambda(il);        
        w_bar=2/3;
        w_bar_last=0;
        cpt=0;
        while (w_bar-w_bar_last)/w_bar>10^-6
            l_bar=(1-w_bar/2)/lambda;
            G=fGoldsteinFactor( l_bar,B,Vr );
            I1=2*trapz(Vr,G.*Vr);
            %I2=2*trapz(Vr,G*l_bar.*Vr.^2./(l_bar^2+Vr.^2));
            I3=2*trapz(Vr,G.*Vr.^3./(l_bar^2+Vr.^2));
            w_bar_last=w_bar;
            w_bar=2/(3*I3)*(I1+I3-sqrt(I1^2-I1*I3+I3^2));
            cpt=cpt+1;
        end
        fprintf('B=%03d - Lambda=%04.1f - Converged after %02d iterations\n',B,lambda,cpt);
        w=w_bar;
        Vlbar(iB,il)=l_bar;
        Vwbar(iB,il)=w_bar;
        Cp(iB,il)=2*w*(1-w/2)*(I1-w/2*I3);
        CT(iB,il)=2*w*(I1-w/2*I3);
    end 
end







% save('OptimalPowerCoeffOkulov.mat')

%%
% load('OkulovCp')
load([PATH.DATA_OUT '/WTTheory/MomentumTheoryActuatorDisk.mat']);
%load('OptimalPowerCoeffOkulov.mat')
% load('data/OptimalPowerCoeffOkulovIt.mat')
% setFigurePath({'C:/work/reports/figs/','C:/work/reports/figsdump/'})
% setMatFigurePath({'C:/work/reports/matfig/','./matfig/'})


close all
figure(1); figure(2);
set(1,'DefaultAxesColorOrder',[0 0 0;],'DefaultAxesLineStyleOrder',':|-.|--|-')
set(2,'DefaultAxesColorOrder',[0 0 0;],'DefaultAxesLineStyleOrder',':|-.|--|-')
legds={};
for iB=1:length(VB)
    legds{end+1}=sprintf('B=%d',VB(iB));
    figure(1)
    hold all
    plot([0 Vlambda],[0 Cp(iB,:)])
    
    figure(2)
    hold all
    plot([0 Vlambda],[0 CT(iB,:)])
    
end
    
figure(1)
box on
ylim([0 0.62])
xlim([0 max(Vlambda)])
plot(vlambdaAD,vCPAD,'k-','LineWidth',2.2)
plot(vlambdaAD,vlambdaAD*0+16/27,'k--','LineWidth',2.2)
%plot(OkulovCp(:,1)-0.02,OkulovCp(:,2),'+')
xlabel('$\lambda$ [.]')
ylabel('$C_P$ [.]')
title('CplambdaOkulov')
legds{end+1}='Momentum theory';
legds{end+1}='Betz';
legend(legds,4)

figure(2)
box on
ylim([0 0.94])
xlim([0 max(Vlambda)])
plot(vlambdaAD,vCTAD,'k-','LineWidth',2.2)
plot(vlambdaAD,vCTAD*0+8/9,'k--','LineWidth',2.2)
xlabel('$\lambda$ [.]')
ylabel('$C_T$ [.]')
title('CTlambdaOkulov')
legend(legds,4)








toc();

M=[ Vlambda' CT(i3,:)' Cp(i3,:)' 1./(Vlbar(i3,:))' Vwbar(i3,:)' Vlbar(i3,:)'];
matrix2latex(M,0,'columnlabels', {'$\lambda$','$C_T$','$C_P$','$1/\bar{l}$', '$\bar{w}$','$\bar{l}$'},'Format',{'%d','%.2f','%.2f','%.1f','%.2f','%.3f'})




%% 

figure(3)
hold all
plot([0 1./Vlbar(i3,:)],[0 1./Vlbar(i3,:)]*0+16/27,'k--','LineWidth',2.2)
plot([0 1./Vlbar(i3,:)],[0 Cp(i3,:)],'k')
plot([0 1./vl_bar],[0 OkulovCp(i3,:)],'k--')
box on
ylim([0 0.62])
xlim([0 20])
xlabel('$1/\\bar{l}$ [.]')
ylabel('$C_P$ [.]')
title('CplbarOkulovB3')


%%

