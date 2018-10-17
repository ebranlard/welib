%% Initialization for the BEM code
%InitClear;
InitDefault;
path(path,'./f_format')
path(path,'./f_aeroelastic')
path(path,'./f_optimal')
path(path,'../CFDPostPro/')

sim='B49'; %<------------- PARAMETER
% sim='NB4_AB00_8';
InitCFDPostPro
%Loading
Algo.Ngrid=80;
fInit('xblade',{XBladeParamFile,PcFile},Opts);
% 
% Rotor.twist=Rotor.twist*0+5.1;
% Rotor.chord=Rotor.chord*0+1.5;

%% calling the BEM code
% Some overwritting of defaultf options
Aero.Wind.V0=[0;0;Vinf];
Rotor.Omega=Omega;
Nacelle.tilt=0;
% Rotor.nB=2;

Algo.NumSect=1;
Simulation.CombinedCase.n=1;
Simulation.CombinedCase.Cases=[Vinf Rotor.Omega Controller.pitch];

%%
% Algo.ReInterp=0;
% Rotor.ProfileSet(1,:)=1;
% Algo.Cl2piAlpha=1; % !!!!!!!!!!!!!!!!!!!!!!!!!!
BEMSimulation


%%

Vrb=Rotor.r;
%%
figure
plot(Vrb,BEM.alpha)
figure
plot(Vrb,BEM.Gamma)

save(sprintf('../VortexLattice/data/%s_BEM.mat',sim),'BEM','Vrb')
% save(sprintf('../VortexLattice/data/%s_BEM2piChord.mat',sim),'BEM','Vrb')
% save(sprintf('../VortexLattice/data/%s_BEM2piChordTwist.mat',sim),'BEM','Vrb')
% save(sprintf('../VortexLattice/data/%s_BEM2pi.mat',sim),'BEM','Vrb')
%save(sprintf('../VortexLattice/data/%s_BEM.mat',sim),'BEM','Vrb')