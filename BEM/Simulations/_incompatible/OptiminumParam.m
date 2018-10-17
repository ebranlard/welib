%% Question 1
InitClear
InitDefault
path(path,'f_optimal')

Algo.ClOverCd=100;
Algo.Cl2piAlpha=1;

%% rotor design
Rotor.r=1:0.1:10;
Rotor.nB=3;
Rotor.rhub=min(Rotor.r);
Rotor.R=max(Rotor.r);
Rotor.ne=length(Rotor.r);
lambda=7;
alpha_d=6;  % alpha(17)=6 deg

Cl=2*pi*alpha_d*pi/180;
Cd=Cl/Algo.ClOverCd;

[a aprime phi c]=getOptimizedParametersClCd(Cl,Cd,lambda,Rotor.nB);
[a2 aprime2 phi2 c2]=getOptimizedParametersClCd(Cl,Cd,lambda-1,Rotor.nB);
[a3 aprime3 phi3 c3]=getOptimizedParametersClCd(Cl,Cd,lambda+1,Rotor.nB);
% phi=(phi+phi2+phi3)/3;
% c=(c+c2+3)/3;


twist=phi-alpha_d;
Rotor.twist=twist';
Rotor.chord=c';

r=Rotor.r;
R=Rotor.R;

%%
rhub=Rotor.rhub;
c0=max(c);
rbar=(r(:)-rhub)/(R-rhub);

options=optimset('Display','iter');
FittedParams=fminsearch(@fFitChord,[0.5 0.5 1 0.05 0.5],options,rbar,c(:)/max(c));

[sse fit]=fFitChord(FittedParams,rbar,c(:)/max(c));
FittedParams
% To check the fit
figure(11)
hold all
plot(rbar,c(:),'*')
plot(rbar,fit*2.2,'r')
%%
options=optimset('Display','iter');
FittedParams=fminsearch(@fFitTwist,[0.5 0.5 1],options,rbar,twist);
[sse fit]=fFitTwist(FittedParams,rbar,c(:)/max(c));
FittedParams

figure(22)
hold all
plot(rbar,twist,'*')
plot(rbar,fit,'r')

%% plotting
figure(14)
hold on
plot([0.1 1], [0.3333 0.3333],'r-.')
plot(r/R,a,'r');
plot(r/R,aprime);

xlabel('r/R')
grid()
box()

legend('Axial induction factor without wake rotation (a=1/3)','Axial induction factor (a)','Tangential induction factor (a`)')
ylabel('Induction factors')
%%
Algo.Ngrid=36;
rfull = linspace(Rotor.rhub,Rotor.R,Algo.Ngrid+1)'; % could be change for more fancy cos
r_mid = (rfull(1:end-1)+rfull(2:end))/2; 
Rotor.r=r_mid;
Rotor.rfull=sort(unique([Rotor.r(:)' Rotor.rhub Rotor.R]));
Rotor.chord         = interp1(r,Rotor.chord,Rotor.r) ;
Rotor.twist         = interp1(r,Rotor.twist,Rotor.r);
Rotor.ne=length(Rotor.r);
Rotor.dr=Rotor.r*0;
Rotor.dr(1)=2*(Rotor.r(1)-Rotor.rhub);
for i=2:length(Rotor.r)
    Rotor.dr(i)=2*(Rotor.r(i)-Rotor.r(i-1)-Rotor.dr(i-1)*0.5);
end


%%
Algo.ClOverCd=100;
Algo.Cl2piAlpha=1;
Algo.TipLoss=1;
Algo.HubLoss=1;
Algo.nbIt=500;
% Algo.Correction='SperaCT';

omega=lambda*6/Rotor.R;

vWS=4:1:12;
vlambda=omega*(Rotor.R)./vWS;

Rotor.SweptArea=pi*(Rotor.R^2-Rotor.rhub^2)
Rotor.twist=Rotor.twist;
Rotor.chord=Rotor.chord;

Cases=zeros(length(vWS),3);
Cases(:,1)=vWS;
Cases(:,2)=Cases(:,2)*0+omega;
Cases(:,3)=Cases(:,3);

Simulation.CombinedCase.n=length(vWS);
Simulation.CombinedCase.Cases=Cases;

Algo.TipLoss=1;
Algo.TipLossMethod='Glauert';
BEMSimulation
BEMGl=BEM;

r=Rotor.r;

%%
figure(3)
clf
hold all
for i=1:length(vWS)
    plot(r/Rotor.R,vBEM(i).Cl)
    plot(r/Rotor.R,vBEM(i).Cd)    
end

figure(2)
clf
hold all
for i=1:length(vWS)
    plot(r/Rotor.R,vBEM(i).a)
end
figure(1)
clf
hold all
for i=1:length(vWS)
    plot(r/Rotor.R,vBEM(i).alpha)
end


figure(4)
clf
hold all
for i=1:length(vWS)
    plot(r/Rotor.R,vBEM(i).Pn)
end