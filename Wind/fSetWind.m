function [ Wind ] = fSetWind(Wind, Sim )
V=Sim.Run.WS;

Wind.V0=[0; 0; V]; % wind speed at hub height
Wind.fV0=@(x)[0; 0; V] ; % wind speed at hub height
end

