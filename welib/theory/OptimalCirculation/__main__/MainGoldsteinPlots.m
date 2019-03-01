%%
InitClear
close all
setFigureLatex(1)
require('OPTIMCIRC','v00')

%%
tic();

%% Comparison with Tibery and Wrench - plots for one lbar and different blades
close all
N=100;  % Number of helical vortex per blade
vx=(1/N:1/N:1);   % [.] 
VB=[2 3 4];
%Vh=[1 1/2 1/3 1/5 1/8 1/12 1/24 1/25 1/36];
Vlbar=[1 1/2 1/4 1/8];
Vlinv=1./Vlbar;

Sty={'r:','b--','g-'};

R=1;
GammaGoldstein=zeros(length(VB),N);
for il=1:length(Vlbar)
    l_bar=Vlbar(il);
    h=2*pi*R*l_bar;
    figure    
    hold all
    grid on, box on
    title(sprintf('h=1/%d',Vlinv(il)));    
    [ KB ] = fCirculationBetz(l_bar,vx);
    plot(vx,KB,'k-','LineWidth',2);
    
    
    text('Interpreter','latex','String',sprintf('$$\\overline{l}=1/%d$$',Vlinv(il)),'Position',[.03 .95],'FontSize',11)
    ylim([0 1])
    xlim([0 1])
    xlabel('$r/R$ [.]')
    ylabel('Normalized Circulation $B\Gamma/hw$ [.]')
    % COmputation of Goldstein function
    for iB=1:length(VB)
        B=VB(iB);
        GammaGoldstein(iB,:) = fGoldsteinFactor( l_bar,B,vx);
        plot(vx/R,GammaGoldstein(iB,:),Sty{iB},'Color',fColrs(iB));
        [r k] = fTiberyWrench( l_bar , B );
        plot(r,k,'ko','MarkerSize',5,'MarkerFaceColor',fColrs(iB),'LineWidth',1)
    end  
end
%%    


%% Comparison Betz Goldstein Prandtl for different lambda
Vlambda=[8 5 2];
B=3;
vx=linspace(0,1,100);
lambda_r=Vlambda(1)*vx;
KB=fCirculationBetz( 1/Vlambda(1),vx);
figure
box on , grid on
hold all
plot(lambda_r,KB,'k-','LineWidth',2.2)
ylim([0 1])
cols=repmat(linspace(0,0.4,3)',[1 3]);
cols=fColrs();
for il=1:length(Vlambda)
    lambda=Vlambda(il);
    lambda_r=lambda*vx;
    l_bar=1/lambda;
    G = fGoldsteinFactor( l_bar,B,vx);
    F = fCirculationPrandtl( lambda,B,vx);
    
    plot(lambda_r,G,'-','Color',cols(il,:))
    plot(lambda_r,F,'--','Color',cols(il,:))    
end
xlabel('$\lambda_r$ [.]')
ylabel('Normalized Circulation $C_\Gamma$ [.]')
legend('Betz','Goldstein','Prandtl','Location','South')
title('Comparison Betz Prandtl Goldstein')

%% Comparison of tip-loss factor for Betz-Goldstein and Prandtl for B=3 and different lambda
Vlambda=[8 5 2];
B=3;
w=1;
vx=linspace(0,1,100);
lambda_r=Vlambda(1)*vx;
figure
box on , grid on
hold all
ylim([0 1])
xlim([0.3 1])
cols=repmat(linspace(0,0.4,3)',[1 3]);
cols=fColrs();
for il=1:length(Vlambda)
    lambda=Vlambda(il);
    lambda_r=lambda*vx/R;
    l_bar=1/lambda;
    G=fTipLossGoldsteinOkulov( l_bar,B,vx);
    F=fTipLossPrandtl( lambda,B,vx);
    
    plot(lambda_r/lambda,G,'-','Color',cols(il,:))
    plot(lambda_r/lambda,F,'--','Color',cols(il,:))    
end
xlabel('$r/R$ [.]')
ylabel('Tip-loss factor $F$ [.]')
legend('Goldstein','Prandtl','Location','South')
title('Comparison Betz Prandtl Goldstein F')


%% four blades Theodorsen p129
vl_inv=[1.4 2 3 4 5];
vl_bar=1./vl_inv;
vx=linspace(0,1,100);
nB=4;
G={};
for il=1:length(vl_bar)
    G{il} =fGoldsteinFactor( vl_bar(il),nB,vx );
end
%%
figure, hold all,grid on ,box on
for il=1:length(vl_bar)
    plot(vx,G{il},'k')
    text('Interpreter','latex','String',sprintf('$$\\overline{l}=1/%s$$',num2str(vl_inv(il))),'Position',[0.6 G{il}(60)+0.035],'FontSize',11)
end
axis equal
ylim([0 1])
xlim([0 1])
xlabel('$x$ [.]')
xlabel('Goldstein factor $K(x)$ [.]')
title('GoldsteinKFourBlades')

%% Two blades  %Theodorsen p128
nB=2;
vl_inv=[0.4 1 2 3 4 5 7 10];
vl_bar=1./vl_inv;
vx=linspace(0,1,100);
G={};
for il=1:length(vl_bar)
    G{il} =fGoldsteinFactor( vl_bar(il),nB,vx );
end
%%
figure, hold all,grid on ,box on
for il=1:length(vl_bar)
    plot(vx,G{il},'k')
    text('Interpreter','latex','String',sprintf('$$\\overline{l}=1/%s$$',num2str(vl_inv(il))),'Position',[0.62 G{il}(62)+0.02],'FontSize',11)
end
axis equal
ylim([0 1])
xlim([0 1])
xlabel('$x$ [.]')
ylabel('Goldstein factor $K(x)$ [.]')
title('GoldsteinKTwoBlades')


%% Influence of N
% Vn=10:10:200;
%  A=zeros(length(Vn),5);
%  for i=1:length(Vn); 
%      A(i,:)=fGoldsteinOkulovF(1/5,3, linspace(0.2,0.95,5)  ,Vn(i)); 
%  end
% 
%  close all
%  plot(Vn,A)
%  


%% Polyfit B=3
R=1;
Vlambda=[16:-2:2];

B=3;
w=1;
vx=linspace(0,R,100);
lambda_r=Vlambda(1)*vx/R;
KB=fCirculationBetz( R/Vlambda(1),vx);
figure
box on , grid on
hold all
% plot(lambda_r,KB,'k-','LineWidth',2)
ylim([0 1])
cols=repmat(linspace(0,0.4,3)',[1 3]);
cols=fColrs();

n_poly=8;
P=zeros(length(Vlambda),n_poly+1);
vsLambda=[];

for il=1:length(Vlambda)
    lambda=Vlambda(il);
    vsLambda{il}=sprintf('$\\mathbf{1/%d}$',lambda);
    
    lambda_r=lambda*vx/R;
    l_bar=1/lambda;
    G=fGoldsteinFactor( l_bar,B,vx );
    F=fCirculationPrandtl( lambda,B,vx );
    
    p=polyfit(vx(:),asin(G(:)),n_poly);
    P(il,:)=p;
%     plot(lambda_r,sin(polyval(p,vx)),'--','Color',cols(il,:))
%     plot(lambda_r,G,'-','Color',cols(il,:))
%     plot(lambda_r,F,'--','Color',cols(il,:))    
end
xlabel('\lambda_r [.]')
ylabel('Normalized Circulation K [.]')
legend('Betz','Goldstein','Prandtl','Location','South')
title('Comparison Betz Prandtl Goldstein')
%%
rowLabels=eval(strcat('{',sprintf('''$a_%d$'',',n_poly:-1:1),'''$a_0$''}'));
rowLabels=rowLabels(end:-1:1)

matrix2latex(P(:,end:-1:1)',0,'rowLabels', rowLabels,'columnLabels', vsLambda, 'format', '%-6.1f', 'alignment', 'r')


%% Polyfit B=2
R=1;
Vlambda=[16:-2:2];

B=2;
w=1;
vx=linspace(0,R,100);
lambda_r=Vlambda(1)*vx/R;
KB=fCirculationBetz( R/Vlambda(1),vx);
figure
box on , grid on
hold all
% plot(lambda_r,KB,'k-','LineWidth',2)
ylim([0 1])
cols=repmat(linspace(0,0.4,3)',[1 3]);
cols=fColrs();

n_poly=8;
P=zeros(length(Vlambda),n_poly+1);
vsLambda=[];

for il=1:length(Vlambda)
    lambda=Vlambda(il);
    vsLambda{il}=sprintf('$\\mathbf{1/%d}$',lambda);
    
    lambda_r=lambda*vx/R;
    l_bar=1/lambda;
    G=fGoldsteinFactor( l_bar,B,vx );
    F=fCirculationPrandtl( lambda,B,vx );
       
    p=polyfit(vx(:),asin(G(:)),n_poly);
    P(il,:)=p;
%     plot(lambda_r,sin(polyval(p,vx)),'--','Color',cols(il,:))
%     plot(lambda_r,G,'-','Color',cols(il,:))
%     plot(lambda_r,F,'--','Color',cols(il,:))    
end
xlabel('\lambda_r [.]')
ylabel('Normalized Circulation K [.]')
legend('Betz','Goldstein','Prandtl','Location','South')
title('Comparison Betz Prandtl Goldstein')
%%
rowLabels=eval(strcat('{',sprintf('''$a_%d$'',',n_poly:-1:1),'''$a_0$''}'));
rowLabels=rowLabels(end:-1:1)

matrix2latex(P(:,end:-1:1)',0,'rowLabels', rowLabels,'columnLabels', vsLambda, 'format', '%-6.1f', 'alignment', 'r')



toc();
