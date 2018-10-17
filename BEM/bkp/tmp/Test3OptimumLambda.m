% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS IS FOR ONE TYPICAL PROFILE
% MAKE IT GENERAL FOR ANY CL/CD RATIO
%% Initialization for the BEM code
InitClear
InitRotor
InitSystemProperties
R=Rotor.R;
r=Rotor.r;5

lambda=8;
alpha_d=6;  % alpha(17)=6 deg

rhub=0.1;
l=([rhub;r]+[r;R])/2;l(1)=rhub;l(length(l))=R;
dr=diff(l);

%% BEM code params
BEMParam.nbIt=500;
BEMParam.alphaCrit=0.01;
BEMParam.relaxation=1;
BEMParam.correction='none';
BEMParam.TipLoss=0;


%% Loop on Number of blades and lambdas
Vlambda=[0.1:0.5:1 1:11];
VB=[1 3 100];
CPlambdaFakeBEM=zeros(length(VB),length(Vlambda));
CPlambdaFakeBEMNoDrag=zeros(length(VB),length(Vlambda));
CPlambda_theory=zeros(length(VB),length(Vlambda));

CdBackup=Rotor.Profiles.Cd;
    
V0=8;
for k=1:length(VB)
    B=VB(k);
    for j=1:length(Vlambda)    
        lambda=Vlambda(j);
        lambda_r=r/R*lambda;   % this is actually lambda_r
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS IS FOR ONE TYPICAL PROFILE
% MAKE IT ALSO GENERAL FOR ANY CL/CD RATIO
        [a aprime phi c beta]=getOptimizedParameters(alpha_d,lambda,B);
        %
        dlambda=dr*lambda/R;
        CPlambda_theory(k,j)=8/lambda^2*sum(aprime.*(1-a).*lambda_r.^3.*dlambda);
%        CPlambda_theory(k,j)=8/lambda^2*trapz(lambda_r,aprime.*(1-a).*lambda_r.^3);
        %fake BEM without DRAG
        Rotor.Profiles.Cd=CdBackup*0;
        BEM = fakeBEM(r,c,beta,a,aprime,B,R,lambda,V0,1.225,phi);
        CPlambdaFakeBEMNoDrag(k,j)=BEM.CP;
        %fake BEM with DRAG
        Rotor.Profiles.Cd=CdBackup;
        BEM = fakeBEM(r,c,beta,a,aprime,B,R,lambda,V0,1.225,phi);
        CPlambdaFakeBEM(k,j)=BEM.CP;
    end

end


%%

figure(31)
hold on
plot(Vlambda,Vlambda./Vlambda*0.5926 ,'k-.')
plot(Vlambda,CPlambda_theory(1,:),'k')
plot(Vlambda,CPlambdaFakeBEMNoDrag(1,:),'b')
plot(Vlambda,CPlambdaFakeBEMNoDrag(2,:),'r')
plot(Vlambda,CPlambdaFakeBEMNoDrag(3,:),'g')

% plotting martin's value on top
% load('B100.mat');
% load 'B3.mat'
% load 'B1.mat'
% plot(B1(:,1),B1(:,2),'b--')
% plot(B3(:,1),B3(:,2),'r--')
% plot(IdealRotorMartin(:,1),IdealRotorMartin(:,2),'g--')

ylim([0.1 0.6])
xlim([0 11])
grid
box
legend('Betz limit','Ideal rotor with wake rotation theory','Optimized rotor : B=1','Optimized rotor : B=3','Optimized rotor : B=100',0)
xlabel('lambda [.]')
ylabel('Cp [.]')
%%
figure(32)
hold on
plot(Vlambda,CPlambdaFakeBEMNoDrag(1,:),'b')
plot(Vlambda,CPlambdaFakeBEMNoDrag(2,:),'r')
plot(Vlambda,CPlambdaFakeBEMNoDrag(3,:),'g')
plot(Vlambda,CPlambdaFakeBEM(1,:),'b--')
plot(Vlambda,CPlambdaFakeBEM(2,:),'r--')
plot(Vlambda,CPlambdaFakeBEM(3,:),'g--')
ylim([0 0.6])
xlim([0 11])
grid
box
legend('Optimized rotor : B=1','Optimized rotor : B=3','Optimized rotor : B=100','Drag taken into account',0)
xlabel('lambda [.]')
ylabel('Cp [.]')

% figure(3)
% hold on
% plot(Vlambda,CPlambda_theory(1,:),'b')
% plot(Vlambda,CPlambdaBEMNoDrag(1,:),'r')
% plot(Vlambda,CPlambdaBEMNoDrag(2,:),'g')
% plot(Vlambda,CPlambdaBEMNoDrag(3,:),'m')
% plot(Vlambda,CPlambdaBEM(1,:),'r:')
% plot(Vlambda,CPlambdaBEM(2,:),'g:')
% plot(Vlambda,CPlambdaBEM(3,:),'m:')
% 
% ylim([0 0.6])
% grid
% box
% legend('Betz theory','B=1','B=3','B=100','With drag')
% xlabel('lambda [.]')
% ylabel('Cp [.]')

%plot(Vlambda,CPlambda)
%plot(Vlambda,CQlambda)