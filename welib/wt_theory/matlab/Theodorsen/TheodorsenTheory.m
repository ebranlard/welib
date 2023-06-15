%%
InitClear
require('OPTIMCIRC','v01');  % v00 is slower but maybe better for small values of lbar
% setFigurePath({'./figs/' , '/work/publications/articles/2012-tiploss-theoretical/figs/'});
setFigurePath({'./figs/' , '/work/publications/phdthesis/figsdump/'});
setFigureTitle(0);
setMatFigure(1);
setFigureLatex(0);
%%
%% verification of the Goldsetin's function - theodorsen page 18 figure 5  - See script Main Goldstein for plot for different l_bar
Rw=10;
l_bar=1/2;
l=l_bar*Rw;
B=2;

nr=100;
vx=linspace(0,1,nr);
[ G ] = fGoldsteinFactor( l_bar,B,vx );
xi=2*trapz(vx,G.*vx);
figure
plot(vx,G)
ylim([0 0.75])
xlabel('r/R [.]')
ylabel('Goldstein Factor K')
title('TheodorsenGoldsteinFactor')  




%% Verification of the mass coefficient - Theodorsen p25 Fig 7 
Rw=10;
l_bar=1/2;
l=l_bar*Rw;
B=2;
nr=20;
nl_bar=20;
vx=linspace(0,1,nr);
vhweird=linspace(0.1,7,nl_bar);
xi=zeros(1,length(vhweird));
for ilb=1:length(vhweird)
    fprintf('.');
    l_bar=vhweird(ilb)/(pi);
    [ G ] = fGoldsteinFactor( l_bar,B,vx );
    xi(ilb)=2*trapz(vx,G.*vx);
end
figure
plot([0 vhweird],[1 xi])
grid on
xlabel('$(U_0+w)/(nD_w) = h/(2R) = \bar{l}/2$ [.]')
ylabel('Mass coefficient $\kappa$ [.]')
title('TheodorsenMassCoefficientB2')



%% Verification of the mass coefficient for different B- Theodorsen p133 -------------------------------- Report plot
Rw=10;
vB=[2 3 4 6 20];
nr=100;
nl_bar=50;
vx=linspace(0,1,nr);
vhweird=linspace(0.01,5,nl_bar); %if too small it collapses
xi=zeros(length(vB),length(vhweird));
for ib=1:length(vB)
    for ilb=1:length(vhweird)
        fprintf('.');
        l_bar=vhweird(ilb)/(pi);
        [ G ] = fGoldsteinFactor( l_bar,vB(ib),vx );
%         if(ilb==1)
%             keyboard
%         end
        xi(ib,ilb)=2*trapz(vx,G.*vx);
    end
end
%

figure,hold all, grid on,box on,legds={};
for ib=1:length(vB)
%      plot([0 vhweird],[1 xi(ib,:)+(1-xi(ib,1))-0.0005],'-','Color',fColrs(length(vB)-ib+1,length(vB)))
     plot([0 vhweird],[1 xi(ib,:)+(1-xi(ib,1))-0.0005],'-','Color',fColrs(ib));
%      plot([0 vhweird],[1 xi(ib,:)])
    legds{end+1}=sprintf('B=%d',vB(ib));
end
warning('Rememebr I''m a little bit cheating here due to umerical integration offset')
l_bar=vhweird/(pi);
xi_inf=1-l_bar.^2.*log(1+1./l_bar.^2);
plot([0 vhweird],[1 xi_inf],'k','LineWidth',2)
xlim([0 5])
legds{end+1}='B=\infty';
xlabel('(U_0+w)/(nD_w) = h/(2R) = l/2 [.]')
ylabel('Mass coefficient \kappa [.]')
title('TheodorsenMassCoefficient')
legend(legds)
% legend(legds,'Location','NorthEast')


%% Tricky Plot - K for different l_bar at different position - Theodorsen p135 B=2 - 136 B=4
Rw=10;
nB=4;
nl_bar=20;
vx=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
vl_bar=linspace(0.1,1,nl_bar);
G=zeros(nl_bar,length(vx));
for ilb=1:length(vl_bar)
    fprintf('.');
    vxx=linspace(0,1,100);
    GG=fGoldsteinFactor( vl_bar(ilb),nB,vxx);
    figure(12),hold all,grid on
    plot(vxx,GG)
    G(ilb,:)  = interp1(vxx,GG,vx);
end
%%
figure,box on,hold all, grid on
mx=repmat(vx,nl_bar,1);
ml=repmat(vl_bar,length(vx),1)';
Ip=1:6;
Id=7:length(vx);
plot(ml(:,Ip),G(:,Ip),'k');
plot(ml(:,Id),G(:,Id),'k--');
xlim([0 1])
xlabel('l_{bar}' )
ylabel('Godstein factor K(x) [.]')
title('GoldsteinFactorDifferentRadialPos')



%% Ideal efficiency against w_bar -Theodorsen p 140

vek=[0 1/10 1/5 2/5 3/5 1];
vw_bar=linspace(0,0.2,50);
figure,hold all,grid on , box on
for ie=1:length(vek)
    eta(ie,:)=(1+vw_bar*(0.5+vek(ie)))./((1+vw_bar).*(1+vek(ie)*vw_bar));
    plot(vw_bar,eta(ie,:),'k')
end
xlim([0 0.225])
ylim([0.9 1])




%% Energy losses - Theodorsen p33 fig 8/
close all
require('VC_LIB_MAT','v02');
require('VC_LIB_C','v03');
Rw=10;
l_bar=1/2;
l=l_bar*Rw;
nB=2;
w=2.4673; % !! Useless but required
w=1;

nhelix=80;
nl_bar=20;
nr=nhelix-1;
npsi=30;

% radial position of helix and CP
vx_helix=linspace(0,1,nhelix);
vx=diff(vx_helix); %!!! Dimension
vx=vx_helix(1:end-1)+vx-vx(1)/2; %CP in between helixes
rhelix=vx_helix*Rw;
% Goldstein Factor
[ K ] = fGoldsteinFactor( l_bar,nB,vx_helix );
xi=2*trapz(vx_helix,K.*vx_helix);


GammaBound=K/nB*w*2*pi*Rw*l_bar;
GammaTrailed=-diff(GammaBound);

% we will consider the difference of gamma to be at the middle points and then will add the smae derivative at the extremity points
% GammaHelix=[GammaTrailed GammaTrailed(end) ];
GammaHelix=[GammaTrailed(1)/2 GammaTrailed ];
figure,plot(vx_helix,K),figure,plot(vx_helix,GammaHelix)

%pitch of the helixes
h=2*pi*l_bar*Rw;

% that's for theoretical solutions where symetry of the rotor is already assumed
vh=repmat(h,1,nhelix); 
vpsi_helix=zeros(1,nhelix);
N=nB; % trick because this code is nasty and need a version 2

Algo = fInitVFilAlgo();
Algo.Method='okulov';
% Algo.Method='line_c';
if(Algo.Method(end)~='N')
    % The computation is done for each blade, so one needs to provided all info for all blades
    vGammaTrailed=repmat(GammaHelix,1,nB); 
    vh=repmat(vh,1,nB); 
    vpsi_helix=reshape(repmat((0:nB-1)'*2*pi/nB,1,nhelix)',1,nB*nhelix);
    vrhelix=repmat(rhelix,1,nB); 
    N=nB*nhelix;
else
    % that's for theoretical solutions where symetry of the rotor is already assumed
    vGammaTrailed=GammaHelix; 
    vh=vh; 
    vpsi_helix=zeros(1,nhelix);
    vrhelix=rhelix; 
    N=nB; % trick because this code is nasty and need a version 2
end
Algo.Helix.bInf=1;
Algo.Helix.nRev=150;
Algo.ntot=10000;
vpsi_CP=linspace(-pi/nB,pi/nB,npsi); %I only take the sector I need, I know it's symetric. Also if npsi is even, I avoid the singularity, if odd then I'm on it and by symmetry it will be everywhere accounted for
% the simulation
vx_CP=[10^-3 vx vx+1]*Rw; % !!! Dimension
% vx_CP=5
% N=2
vz=(Algo.Helix.nRev/2+0.0000)*h;
[ui Xcp Ycp Zcp Helix]=fUi_HelixN(N,vx_CP,vpsi_CP,vz,vGammaTrailed,vrhelix,vh,vpsi_helix,Algo); % the z position doesnot matter since we are dealing with an infinite helix here
ui
%%
% figure(12),hold all
% for ih=1:length(Helix)
%     plot3(Helix{ih}{1}(:,1), Helix{ih}{1}(:,2), Helix{ih}{1}(:,3)) 
% end

% ui goes length(vx) length(vpsi) 1 [r theta z]
Rmax=max(vx_CP);
  figure
    h3 = polar([0 2*pi], [0 Rmax],'');delete(h3)
    hold on
    contourf(Xcp,Ycp,squeeze(ui(:,:,1,3)),30,'LineStyle','none')
    colorbar
    title('ax')
  figure
    h3 = polar([0 2*pi], [0 Rmax],'');delete(h3)
    hold on
    contourf(Xcp,Ycp,squeeze(ui(:,:,1,2)),30,'LineStyle','none')
    colorbar
    title('tan')
  figure
    h3 = polar([0 2*pi], [0 Rmax],'');delete(h3)
    hold on
    contourf(Xcp,Ycp,squeeze(ui(:,:,1,1)),30,'LineStyle','none')
    colorbar
title('rad')
%%
% Area far wake
Ir=1:sum(vx_CP<Rw);

F=pi*Rw^2 /nB  %!!! By symmetry, since I use symmetry for integration, I divide by B
inte=repmat(vx_CP(Ir),length(vpsi_CP),1)';
F=trapz(vpsi_CP,trapz(vx_CP(Ir),inte))

% quad dlbquad
%
ur2=squeeze(ui(:,:,:,1)).^2;
ut2=squeeze(ui(:,:,:,2)).^2;
uz2=squeeze(ui(:,:,:,3)).^2;


% trying to add some point to make integration better
vr0=[0 vx_CP(Ir) Rw];
uz20=[repmat(0,1,length(vpsi_CP));uz2(Ir,:); repmat(0,1,length(vpsi_CP))];
epsz0=1/(2*F)*trapz(vpsi_CP,trapz(vr0,uz20.*repmat(vr0,length(vpsi_CP),1)'))

epsr=1/(w^2*F)*trapz(vpsi_CP,trapz(vx_CP(Ir),ur2(Ir,:).*inte))
epst=1/(w^2*F)*trapz(vpsi_CP,trapz(vx_CP(Ir),ut2(Ir,:).*inte))
epsz=1/(w^2*F)*trapz(vpsi_CP,trapz(vx_CP(Ir),uz2(Ir,:).*inte))
fprintf('Sum eps: %.4f xi: %.4f\n',epsr+epst+epsz,xi);

depsr=1/(w^2*F)*trapz(vpsi_CP,ur2')/2*Rw^2;  % depsilon/d(x^) = depsilon/dr * R^2 /(2r)   but the r cancels out due to the expression r dr dtheta
depst=1/(w^2*F)*trapz(vpsi_CP,ut2')/2*Rw^2;
depsz=1/(w^2*F)*trapz(vpsi_CP,uz2')/2*Rw^2;


depst2=[depst(1:length(vx)+1) 0 depst(length(vx)+2:end)];

vx_CP2=[vx_CP(1:length(vx)+1) Rw vx_CP(length(vx)+2:end)];
depst2=[depst(1:length(vx)+1) 0 depst(length(vx)+2:end)];
depsz2=[depsz(1:length(vx)+1) 0 depsz(length(vx)+2:end)];

%%
figure,hold all,grid on,box on
plot(vx_CP.^2./Rw^2,depsr,'k-','LineWidth',2.2)
plot(vx_CP.^2./Rw^2,depst,'k:','LineWidth',1.5)
plot(vx_CP.^2./Rw^2,depsz,'k--','LineWidth',1.5)
xlim([0 3])
ylim([0 0.36])
set(gca,'ytick',0:0.04:0.36)
legend('r','t','z')
%% testint
vx=[0 vx 1]
inte=repmat(vx*Rw,length(vpsi_CP),1)';
trapz(vpsi_CP,trapz(vx*Rw,inte))
trapz(vx*Rw,trapz(vpsi_CP,inte'))
%%

a_ll=squeeze(-ui(:,1,:,3))'/V0;



[ ui ] = fUi_HelixNTheories( 'okulov',Gamma,r,r0 ,l,psih,nB, bWT,bHalf);


%% Loss function epsilon/xi for an infinite number of blades - Thodorsen p34
l_bar=linspace(0,2,100);
epsz=1+l_bar.^2./(1+l_bar.^2)-2*l_bar.^2.*log(1+1./l_bar.^2);
epst=-l_bar.^2./(1+l_bar.^2)+l_bar.^2.*log(1+1./l_bar.^2);
xi_inf=1-l_bar.^2.*log(1+1./l_bar.^2);
epsoverxi=epsz./xi_inf;
figure,hold on,box on,grid on
plot(l_bar,epsz,'k')
plot(l_bar,epst,'k')
plot(l_bar,xi_inf,'k')
plot(l_bar,epsoverxi,'k')
legend('\epsilon','\epsilon_t','\kappa','\epsilon/\kappa')


%% Loss function estimate - Theodorsen p37  --------------------------------------Report plot
Rw=10;
l=l_bar*Rw;
nB=4; % B=2 or B4
nB=3; % for report
nr=100;
nl_bar=100;
vx=linspace(0,1,nr);
vl_bar=linspace(0.01/pi,5/pi,nl_bar);  % the new algo does not support small values as well
xi=zeros(1,length(vl_bar));
for ilb=1:length(vl_bar)
    fprintf('.');
    [ G ] = fGoldsteinFactor( vl_bar(ilb),nB,vx );
    xi(ilb)=2*trapz(vx,G.*vx);
end
dxidl=diff(xi)./diff(vl_bar);  % slope for B=2 around 3pi/5 and pi/25 for B=4 at 0
% dxidl=[dxidl(1) dxidl];
dxidl=[0 dxidl];
epsoverxi=1+0.5*vl_bar./xi.*dxidl;
epsz=epsoverxi.*xi;


figure, grid on, box on, hold all
plot([0 vl_bar*pi],[1 xi],'k-','LineWidth',2)
plot([0 vl_bar*pi],[1 epsz],'k-')
plot([0 vl_bar*pi],[1 epsoverxi],'k--')
xlim([0 5])
ylim([0 1])
legend('\kappa','\epsilon','\epsilon/\kappa')
grid on
xlabel('(U_0+w)/(nD_w) = h/(2R) = l\pi [.]')
ylabel('[.]')
title(['TheodorsenLossEstimateB' num2str(nB)])

%%


