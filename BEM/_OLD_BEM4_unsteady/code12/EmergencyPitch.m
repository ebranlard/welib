function theta=EmergencyPitch(t,tpitch,tstop)
    global Controller Algo
    if(t<tstop)
        theta=PitchControl(t);
    else
        if(Controller.lastpitch<35)
            theta=Controller.pitch+35*Algo.dt/tpitch;
        else
            theta=35;
        end
    end
end