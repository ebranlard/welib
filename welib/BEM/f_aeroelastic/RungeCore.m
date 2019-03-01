
%% time loop
tic();
for i=INDEX
%     if(abs(mod(t,2*pi/omega))<=dt/2)
%         disp(['Turn ' num2str(round(t/(2*pi/omega))) ' - t=' num2str(t)])
%     end
    if(mod(i,10)==1)
         disp(['Iteration ' num2str(i) ' - t=' num2str(t(i))])
    end
    Controller.pitch=Controller.fpitch(t(i)); 
    
    if(isfield(Algo,'ConstantOmega')==1)
        x(:,i+1)=[0; x(2,i)+v(2,i)*dt ;0];
        v(:,i+1)=v(:,i);
        t(i+1)=t(i)+dt;
    else
        %%%% Runge Kutta Nystrom scheme
        %%% Estimates at t+dt/2
        A=dt/2*a(:,i);
        b=dt/2*(v(:,i)+0.5*A);
        % speed changed to V+A, position changed to x+b
        B=dt/2.*solveAcceleration( t(i)+dt/2 , x(:,i)+b , v(:,i)+A, NoUpdate);
        % speed changed to V+B
        C=dt/2 *solveAcceleration(  t(i)+dt/2 , x(:,i)+b , v(:,i)+B, NoUpdate);

        %%% Estimates at t+dt
        d=dt*(v(:,i)+C);
        % speed changed to v+2*C, position changed to x+d
        D=dt/2*solveAcceleration(  t(i)+dt , x(:,i)+d , v(:,i)+2*C ,NoUpdate);

        %%% Final results
        t(i+1)=t(i)+dt;
        x(:,i+1)=x(:,i)+dt.*(v(:,i)+1/3.*(A+B+C));
        v(:,i+1)=v(:,i)+1/3*(A+2*B+2*C+D);
    end
    if(Algo.DontRotate)
        v(2,i+1)=0;
        a(2,i+1)=0;
        x(2,i+1)=0;
        v(3,i+1)=0;
        a(3,i+1)=0;
        x(3,i+1)=0;
    end
    % final speed and positions known
    a(:,i+1)=solveAcceleration( t(i+1) , x(:,i+1) , v(:,i+1), Update);
    
    %%%% Controller
    Controller.lastpitch=Controller.pitch;
    if(Controller.pitch==35)
        Algo.DontRotate=1
    end
    
    
    %%%% Time Parameters storage
    Pitch(i) = Controller.pitch;
    Torque(i) = Aero.Torque;
    Thrust(i) = Aero.Thrust;
    Power(i)  = Aero.Power;
    HubV(i)=Aero.Wind.V0(3);
    MG(i)  = Generator.fMoment(v(2,i));    
    Flap(:,i)  = Aero.Flap;   
        
    %%%% Auto Stop when power converges
    if(t(i)>25)
        if(isfield(Algo,'AutoStop')==1 && abs((Power(i)-Power(i-1))/Power(i))<0.001  && abs((Power(i)-Power(i-2))/Power(i))<0.001 )
            disp('Convergence clause reached -> Breaking.')
            break;
        end
    end
    if(isfield(Algo,'breakdown')==1)
        disp('RungeCore caught BEM crash -> Breaking.');
        break;
    end
end
Torque(i+1) = Torque(i);
Pitch(i+1) = Pitch(i);
Thrust(i+1) = Thrust(i);
Power(i+1)  = Power(i);
HubV(i+1)  = HubV(i);
MG(i+1)  = MG(i);
Flap(:,i+1)  = Flap(i);   
toc()