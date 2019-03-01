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
Algo.HubLoss=0;
Algo.TIDrag=1;
Simulation.CombinedCase.n=1;
Simulation.CombinedCase.Cases=[5.5 700.2817*2*pi/60 2];
%%
Algo.Cl2piAlpha=1;
Algo.ClOverCd=100;
Rotor.twist=Rotor.twist*0+5.1;
BEMSimulation




%%
Vrb=Rotor.r;
figure
plot(Vrb,BEM.alpha)
grid on
figure
plot(Vrb,BEM.Gamma)
grid on

% save('../VortexLattice/data/FlatTUD_BEM2pi.mat','BEM','Vrb')
% save('../VortexLattice/data/FlatTUD_BEM.mat','BEM','Vrb')