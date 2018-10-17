function g=solveAcceleration(t,x,v,update)
    global Generator Shaft Nacelle Tower Rotor Aero Controller Algo
    
    %%% Time dependent parameters
    Controller.yaw=Controller.fyaw(t);
    %Controller.pitch=Controller.fpitch(t);    
    Aero.Wind.V0=Aero.Wind.fV0(t);
    Algo.t=t;
    
    %% Matrices
    [M K C]=getMatrices();

    %% Generalized forces
    %Calling the BEM code for this new positions and speeds
    [BEM]=fBEM(x,v,update);
    %Generalized forces
    GF=zeros(3,1);
    GF(1)=BEM.Thrust;
    GF(2)=BEM.Torque-Generator.fMoment(v(2));
    GF(3)=-Generator.fMoment(v(2)+v(3));
    GF(4:6)=BEM.GF(1,:);
    GF(7:9)=BEM.GF(2,:);
    GF(10:12)=BEM.GF(3,:);
   
    if(Algo.Break)
        if(v(2)>0.001)
             GF(2)=-10^3;
        else
            GF(2)=0;
        end
        C(2,2)=10^10;
        if(Power<0)
            Aero.Power=0;
        end
    end
    
    % System solving for the acceleration
    g = inv(M)*(GF-K*x-C*v) ;
end
