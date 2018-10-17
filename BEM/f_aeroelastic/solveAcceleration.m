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
    GF(3)=-Generator.fMoment(v(2));
    
    % 2DOF
    if(Algo.TwoDOF)
        M=M(1:2,1:2);
        K=K(1:2,1:2);
        C=C(1:2,1:2);
        GF=GF(1:2);
    end
    
    % System solving for the acceleration
    g = inv(M)*(GF-K*x-C*v) ;
end
