%% Initialization for the BEM code
InitClear;
InitDefault;
path(path,'./f_format')
path(path,'./f_aeroelastic')
path(path,'./f_optimal')

Algo.Ngrid=36;
fInit('xblade',{'C:/work/data/Aero/FlatTUD/FlatTUD_param.txt','C:/work/data/Aero/FlatTUD/NACA0012.pc'},Opts)

%% calling the BEM code 
Algo.Steady=1;
Algo.NumSect=1;
Algo.nbIt=500;
%Algo.alphaCrit=0.001;
%Algo.relaxation=0.8;
%Algo.correction='Spera';
Algo.YawModel=0;
Algo.TipLoss=1;
Algo.HubLoss=1;
Algo.TIDrag=1;


vlambda=1:14;
omega=700*2*pi/60;
vWS=Rotor.R*omega./vlambda;
Rotor.SweptArea=pi*(Rotor.R^2-Rotor.rhub^2)

Cases=zeros(length(vWS),3);
Cases(:,1)=vWS;
Cases(:,2)=Cases(:,2)+omega;

Simulation.CombinedCase.n=length(vWS);
Simulation.CombinedCase.Cases=Cases;
%%
Algo.Cl2piAlpha=1;
% Rotor.twist=Rotor.twist*0+5.1;
BEMSimulation




%%
Vrb=Rotor.r;
figure
plot(Vrb,BEM.alpha)
grid on
figure
plot(Vrb,BEM.Gamma)
grid on

save('../VortexLattice/data/FlatTUD_BEM2pi.mat','BEM','Vrb')
% save('../VortexLattice/data/FlatTUD_BEM.mat','BEM','Vrb')