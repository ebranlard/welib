%% Initialization for the BEM code
loading;
V0=8; 
R=3;
B=3;
rhub=0.1;
    
%% BEM code parameters
Param.nbIt=5000;
Param.alphaCrit=0.05;
Param.relaxation=1;
Param.correction='none';

%% question 4 first part
R=3;
omega=21.33;
r=2.6;
c=0.1;  
dr=0.1;
Vspeed=5:9;
Vbeta=-10:0.5:20;  
Cpbeta=zeros(length(Vspeed),length(Vbeta));
betaV0=zeros(1,length(Vspeed));
for jj=1:length(Vspeed)
    V0=Vspeed(jj);
    lambda=omega*R/V0;
    x=r/R*lambda;
    for ii=1:length(Vbeta)
        a=0.2;
        aprime=0.01;
        beta=Vbeta(ii);
        [BEM CP CT CQ] =BEMfunction(Param,Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho);
        Cpbeta(jj,ii)=BEM.Cp;
    end
    index=whichmax(Cpbeta(jj,:));
    betaV0(jj)=Vbeta(index);
end
betaopt=sum(betaV0.*wblpdf(Vspeed,6,2))/sum(wblpdf(Vspeed,6,2))


figure(3)
hold on
plot(Vbeta,Cpbeta(1,:),'b')
plot(Vbeta,Cpbeta(2,:),'r')
plot(Vbeta,Cpbeta(3,:),'g')
plot(Vbeta,Cpbeta(4,:),'m')
plot(Vbeta,Cpbeta(5,:),'k')
grid()
box()
ylabel('C_p [.]')
xlabel('Twist angle [deg]')
legend('V=5 m/s','V=6 m/s','V=7 m/s','V=8 m/s','V=9 m/s')


%% Weibull
% 
% plot(0:0.1:25,wblpdf(0:0.1:25,6,2))
% xlabel('Wind speed [m/s]')
% ylabel('Density of probability [.]')
%     grid()
%     box()