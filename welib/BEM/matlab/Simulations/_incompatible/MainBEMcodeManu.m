%% Initialization for the BEM code
InitClear;
InitDefault;
path(path,'./f_format')
path(path,'./f_aeroelastic')
path(path,'./f_optimal')

Algo.Ngrid=45;
fInit('xblade',{'C:/work/data/Aero/Manu/Manu_param.txt','C:/work/data/Aero/Manu/NACA0012.pc'},Opts)
%%
Algo.ClOverCd=100;
Algo.Cl2piAlpha=1;
Algo.TipLoss=0;
Algo.HubLoss=1;
Algo.nbIt=500;

omega=4.2;

vWS=4:1:14;
vlambda=omega*(Rotor.R)./vWS;

Rotor.SweptArea=pi*(Rotor.R^2-Rotor.rhub^2)
Rotor.twist=Rotor.twist;
Rotor.chord=Rotor.chord;


Cases=zeros(length(vWS),3);
Cases(:,1)=vWS;
Cases(:,2)=Cases(:,2)+omega;
Cases(:,3)=Cases(:,3);

Simulation.CombinedCase.n=length(vWS);
Simulation.CombinedCase.Cases=Cases;

Algo.TipLoss=1;
Algo.TipLossMethod='Glauert';
BEMSimulation
BEMGl=BEM;

r=Rotor.r;

%
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


figure(6)
clf
hold all
for i=1:length(vWS)
    plot(r/Rotor.R,vBEM(i).Cl./vBEM(i).Cd)
end
%%
figure(7)
clf
hold all
for i=1:length(vWS)
    plot(r/Rotor.R,vBEM(i).Gamma)
end