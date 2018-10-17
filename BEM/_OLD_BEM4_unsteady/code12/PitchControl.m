function theta=PitchControl(t)
    global Aero Controller Algo
    theta=Controller.pitch;
    if(t>Algo.dt)
        tau1=1.1;
        tau2=Algo.dt*1.1;
        theta_c=0.02*(Aero.Power/1000-2000)/(1+(Controller.lastpitch)/4.6)*tau1+Controller.lastpitch;        
        theta=max(theta_c+(Controller.lastpitch-theta_c)*exp(-Algo.dt/tau2),0);
    end
end