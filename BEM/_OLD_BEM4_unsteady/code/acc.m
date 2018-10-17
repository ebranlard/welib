function g=solveAcceleration(t,x,v,update)
    global Generator Shaft Nacelle Tower Rotor Aero

    %%% Time dependent parameters
%     VelocityParams.V0
    
    %% Matrices
    [M K C]=getMatrices();

    %% Generalized forces
    %Calling the BEM code for this new positions and speeds
    [BEM]=fBEM(x,v,update)
    %Generalized forces
    GF=zeros(3,1);
    GF(1)=BEM.Thrust;
    GF(2)=BEM.Torque-Generator.fMoment(v(2));
    GF(3)=Generator.fMoment(v(2));
    
    % System solving for the acceleration
    g = inv(M)*(GF-K*x-C*v) ;
end
