function fInitWTPerf(Files,Opts)
global Rotor Controller Aero Simulation Algo
fReadWTPerf(Files{1});

Rotor.thickness_rel=zeros(1,Rotor.ne);

Rotor.Omega=Simulation.Parametric.Omega(1);
Controller.pitch=Simulation.Parametric.Pitch(1);
Aero.Wind.V0=[0;0;Simulation.Parametric.WS(1)];
Algo.relaxation=0.3;