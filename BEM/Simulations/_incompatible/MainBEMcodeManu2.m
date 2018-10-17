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

Simulation.CombinedCase.n=1;
Simulation.CombinedCase.Cases=[6 omega 0];
%%
BEMSimulation




%%
Vrb=Rotor.r;
figure
plot(Vrb,BEM.alpha)
grid on
figure
plot(Vrb,BEM.Gamma)
grid on

save('../VortexLattice/data/Manu_BEM2pi.mat','BEM','Vrb')
% save('../VortexLattice/data/FlatTUD_BEM.mat','BEM','Vrb')