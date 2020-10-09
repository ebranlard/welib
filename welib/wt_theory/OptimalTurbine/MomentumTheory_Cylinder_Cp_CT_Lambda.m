%%
InitClear
% setFigurePath('/work/publications/book/figsdump/')
% I explain this a bit in my book section Maximum Power extraction 1D momentum+rotation

setFigurePath('/work/publications/articles/2014-cylinder-superp/figs/')
setMatFigurePath('/work/publications/articles/2014-cylinder-superp/matfigs/')

setFigureFont('15')

%
R=100;
nr=200;

n2=10000; % for second method
n3=100000; % for third method


options = optimset('TolX',1e-11);


%% Vortex cylinder way
% vlambda=unique([linspace(0.1,2,12) linspace(2,10,30)]);
vlambda=linspace(0.1,10,40);
CP=zeros(1,length(vlambda));
CPb=zeros(1,length(vlambda));
CT=zeros(1,length(vlambda));
aprime=zeros(1,length(vlambda));
kmax=zeros(1,length(vlambda));
a=zeros(1,length(vlambda));


for il=1:length(vlambda)
    lambda=vlambda(il);
%     kmax=fzero(@(k) -k^2-3*k*lambda^2+(2+sqrt(4+k*(-4-k/lambda^2)))*lambda^2 ,options); %martin 8.1 combined from 4.32 and 4.38
    kmax=fzero(@(k) -8*lambda^4 + (-3*lambda^2 + 9*lambda^4)*k + 6*lambda^2*k^2 + k^3,0.5,options);
    % --
    CP(il)=kmax*(1+0.5*(-1+sqrt(1-kmax*(1+kmax/(4*lambda^2)))));

    % --- NEW 
    a0=0.5*(1-sqrt(1+lambda^2)*sin(1/3*atan(1/lambda))); % Eq(48) in branlard:2014 cylinder superp
    kmax=2*lambda*(sqrt(lambda^2+4*a0*(1-a0))-lambda );
    CP2(il)=kmax*(1+0.5*(-1+sqrt(1-kmax*(1+kmax/(4*lambda^2)))));
%     CT(il)=kmax*(1+0.5*(-1+sqrt(1-kmax*(1+kmax/(4*lambda^2)))));
end


%% Momentum Theory Way - Deterministic Way
tic()
disp('Third method...')
a3=linspace(0.25+10^-9,1/3,n3); % why do i start at 0.25, because of the expression for aprime
aprime3=(1-3*a3)./(4*a3-1);
vlambda3=sqrt(a3.*(1-a3)./(aprime3.*(1+aprime3)));
[vlambda3,Isort]=sort(vlambda3);
CP3=8./(vlambda3.^2).*cumtrapz(vlambda3,aprime3(Isort).*(1-a3(Isort)).*vlambda3.^3);
CT3=8./(vlambda3.^2).*cumtrapz(vlambda3,a3(Isort).*(1-a3(Isort)).*vlambda3);
a3=a3(Isort);
aprime3=aprime3(Isort);
toc()




%% CP
InitCell
N4=load('OkulovSorensen2010_N4');
N3=load('OkulovSorensen2010_N3');
N1=load('OkulovSorensen2010_N1');


figure,hold on,grid,box
plot(vlambda,vlambda*0+16/27 ,'k-','LineWidth',2)
plot(vlambda,CP,'d','Color',fColrs(1))
plot(vlambda3,CP3,':','Color',fColrs(2),'LineWidth',3)
plot(N4.OkulovSorensen2010_N4(:,1),N4.OkulovSorensen2010_N4(:,2),'--','Color',fColrs(4))
% plot(N3.OkulovSorensen2010_N3(:,1),N3.OkulovSorensen2010_N3(:,2),'k')
% plot(N1.OkulovSorensen2010_N1(:,1),N1.OkulovSorensen2010_N1(:,2),'k')
set(gca(),'YTick',[0,0.2,0.4 16/27])
set(gca(),'YTickLabel',{'0','0.2','0.4','16/27'})
legend('Betz-Joukowski limit','Current (cylindrical)','Stream-tube Theory (Glauert)','4 Helices, h\neqEq. 13, (Okulov et al.)','Location','SouthEast')
xlim([0 10])
ylim([0 0.63])
xlabel('Tip speed ratio \lambda [.]')
ylabel('Power Coefficient C_P [.]')
title('OptimalCPWithVortexCylinder')



%% CT
% figure,hold on,grid,box
% plot(vlambda,vlambda*0+8/9 ,'k--')
% plot(vlambda3,CT3,'r--')
% plot(vlambda,CT,'b--')
% legend('Betz limit','Classical Momentum Theory','Vortex Cylinder',0)
% xlim([0 15])
% xlabel('Tip speed ratio \lambda [.]')
% ylabel('C_T [.]')
% title('MomentumTheoryActuatorDiskOptimalCT')

% %% a
% figure,hold on,grid,box
% plot(vlambda,vlambda*0+1/3 ,'k--')
% plot(vlambda3,a3,'r--')
% legend('Betz limit','Classical Momentum Theory','Vortex Cylinder',0)
% xlim([0 15])
% ylim([0 0.35])
% xlabel('Tip speed ratio \lambda [.]')
% ylabel('a [.]')
% title('MomentumTheoryActuatorDiskOptimala')
% 
% %% aprime
% figure,hold on,grid,box
% plot(vlambda3,aprime3,'r--')
% ylim([0 0.1])
% legend('Betz limit','Classical Momentum Theory','Vortex Cylinder',0)
% xlim([0 15])
% title('MomentumTheoryActuatorDiskOptimalaprime')
% xlabel('Tip speed ratio \lambda [.]')
% ylabel('a'' [.]')
% 


save('MomentumTheoryCylinder_CP.mat','vlambda','CP','vlambda3','CP3')


