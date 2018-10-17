%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TWIST OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization for the BEM code
InitClear
InitRotor
InitSystemProperties
R=Rotor.R;
r=Rotor.r;
omega=Rotor.Omega;

%% deciding a chord
%deciding a lambda
Vspeed=5:9;
Vlambda=omega*R./Vspeed;
lambdachosen=sum(Vlambda .*wblpdf(Vspeed,6,2))/sum(wblpdf(Vspeed,6,2))
%% Determining the optimum chord, for this lambda and using it
lambda=lambdachosen;
alpha_d=6;  % alpha(17)=6 deg
[a aprime phi c beta]=getOptimizedParameters(alpha_d,lambda,Rotor.nB);
chord=c;


%% BEM code parameters
Param.nbIt=5000;
Param.alphaCrit=0.05;
Param.relaxation=0.3;
Param.correction='Spera';

%%
Vspeed=5:9;
Vbeta=-10:0.5:30;  
betaV0=zeros(length(r),length(Vspeed)); %optimized twist for each position and wind speed
betaopt=zeros(1,length(r));             %optimized twist for each position

Cpbeta=zeros(length(Vspeed),length(Vbeta));
    for jj=1:length(Vspeed) %loop different wind speed
        V0=Vspeed(jj);
        lambda=omega*R/V0;
        x=r/R*lambda;
        for ii=1:length(Vbeta) %loop different beta
            a=0.2;        %initialization for BEM
            aprime=0.01;  %initialization for BEM
            beta=Vbeta(ii);
            %BEM calculation
            %[BEM CP CT CQ] = BEMfunction(Param,Profile,r,dr,x,c,beta,a,aprime,B,R,lambda,V0,rho); 
            BEM = fBEMsteady(0,0,0);
            Cpbeta(jj,ii)=BEM.Cp;
        end %loop different beta
        index=whichmax(Cpbeta(jj,:));
        betaV0(kk,jj)=Vbeta(index);
end  %loop different wind speed
betaopt(kk)=sum(betaV0(kk,:).*wblpdf(Vspeed,6,2))/sum(wblpdf(Vspeed,6,2));

betaopt
%% comparison with optimized value at the lambda chosen
lambda=lambdachosen;
alpha_d=6;  % alpha(17)=6 deg
r=rsection;
x=r/R*lambda;   % this is actually lambda_r
[a aprime phi c beta ]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x);

figure(41)
hold on
plot(r/R,betaopt,'b')
plot(r/R,beta,'r')
xlabel('r/R [.]')
grid()
box()
ylabel('Twist angle [deg]')
legend('Twist - optimum','Twist - theory ')
