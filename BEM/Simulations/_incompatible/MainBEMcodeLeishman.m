%% Initialization for the BEM code
InitClear;
InitDefault;
path(path,'./f_format')
path(path,'./f_aeroelastic')
path(path,'./f_optimal')

Algo.Ngrid=36;
fInit('xblade',{'C:/work/data/Aero/Leishman/Leishman_param.txt','C:/work/data/Aero/Leishman/NACA0012.pc'},Opts)

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
omega=20/3;
vWS=Rotor.R*omega./vlambda;


Rotor.SweptArea=pi*(Rotor.R^2-Rotor.rhub^2)
Rotor.twist=Rotor.twist-1;

% vWS=[2 2.5 4 5 8 10 ];
% vlambda=omega*Rotor.R./vWS;

Cases=zeros(length(vWS),3);
Cases(:,1)=vWS;
Cases(:,2)=Cases(:,2)+omega;

Simulation.CombinedCase.n=length(vWS);
Simulation.CombinedCase.Cases=Cases;
%%
 Algo.Cl2piAlpha=1;
BEMSimulation


%%
figure
plot([0 vlambda],[0 CP])
figure
plot([0 vlambda],[0 CT])

