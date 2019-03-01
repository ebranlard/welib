%% Initialization for the BEM code
loading;
V0=8; 
R=3;
rhub=0.1;
alpha_d=6;  % alpha(17)=6 deg
r= [0.3 0.5 1.0 1.5 2.0 2.5 2.6 2.7 2.8 2.9 2.95 ];
%length of elements
l=([rhub r]+[r R])/2;l(1)=rhub;l(length(l))=R;
dr=diff(l);

% BEM code parameters
Param.nbIt=5000;
Param.alphaCrit=0.05;
Param.relaxation=1;
Param.correction='Glauert';

Vlambda=[0.001:0.5:4.001 5:2:20 20]
VB=[1 3 100];
CPlambdaFakeFEM=zeros(length(VB),length(Vlambda));
CPlambdaFakeFEMNoDrag=zeros(length(VB),length(Vlambda));
CPlambdaFEM=zeros(length(VB),length(Vlambda));
CPlambdaFEMNoDrag=zeros(length(VB),length(Vlambda));
CPlambda_theory=zeros(length(VB),length(Vlambda));
for k=1:length(VB)
    B=VB(k);
    for j=1:length(Vlambda)
        lambda=Vlambda(j);
        x=r/R*lambda;   % this is actually lambda_r
        [a aprime phi c beta]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x);
        %% definition of elements length
        dlambda=dr*lambda/R;
        CPlambda_theory(k,j)=8/lambda^2*sum(aprime.*(1-a).*x.^3.*dlambda);
        %fake BEM
        
        [BEM CP CT CQ] = fakeBEM(Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho,phi);   
%         [BEM CP CT CQ] =BEMfunction(Param,Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho)
        CPlambdaFakeBEM(k,j)=CP;
        [BEM CP CT CQ] = fakeBEM(ProfileNoDrag,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho,phi);   
        CPlambdaFakeBEMNoDrag(k,j)=CP;
        %BEM        
%         [BEM CP CT CQ] = BEMfunction(Param,Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho);   
%         CPlambdaBEM(k,j)=CP;
%         [BEM CP CT CQ] = BEMfunction(Param,ProfileNoDrag,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho);   
%         CPlambdaBEMNoDrag(k,j)=CP;
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
ylim([0 0.6])
grid
box
legend('Betz limit','Ideal rotor with wake rotation theory','Optimized rotor : B=1','Optimized rotor : B=3','Optimized rotor : B=100')
xlabel('lambda [.]')
ylabel('Cp [.]')
%%
figure(32)
hold on
%plot(Vlambda,CPlambda_theory(1,:),'k')
plot(Vlambda,CPlambdaFakeBEMNoDrag(1,:),'b')
plot(Vlambda,CPlambdaFakeBEMNoDrag(2,:),'r')
plot(Vlambda,CPlambdaFakeBEMNoDrag(3,:),'g')
plot(Vlambda,CPlambdaFakeBEM(1,:),'b--')
plot(Vlambda,CPlambdaFakeBEM(2,:),'r--')
plot(Vlambda,CPlambdaFakeBEM(3,:),'g--')
%ylim([0 0.6])
grid
box
legend('Optimized rotor : B=1','Optimized rotor : B=3','Optimized rotor : B=100','Drag taken into account')
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