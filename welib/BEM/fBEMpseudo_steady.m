function BEM=fBEMpseudo_steady()
global Shaft
global Nacelle
global Tower
global Aero
global Environment
global Generator
global Controller
global Algo
global Rotor
global Profiles
global Simulation
global Spec


Algo.dt=0.1;
t_max=50;
Vtime=0:Algo.dt:t_max;

%% Initializations
nt=length(Vtime);
Omega=Rotor.Omega;
nB=Rotor.nB;
ne=Rotor.ne;
dt=Algo.dt;

psi=0; %[rad] !!!
Power=zeros(1,nt);

%% Time loop
Plast=-1;
for idT=1:length(Vtime)
    t=Vtime(idT);
    %%% Time dependent parameters
    x=[0 psi];    %psi [rad]
    v=[0 Omega];  %Omega [rad/s]
    BEM=fBEM(x,v,1);
 
    %%% Rotation of the blades
    psi=mod(psi+Omega*dt,2*pi) ; %[rad]
    
    %%%% Time Parameters storage
    Power(idT)  = Aero.Power;
    if(idT>10)
        Pmean=mean(Power((idT-9):idT));
        if(abs((Plast-Pmean)/Plast)<Algo.swTol)
            disp(['Converged after : ',num2str(idT),' iterations'])
            break;
        end
        Plast=Pmean;
    end
    % Convergence?
end % time loop
if(idT==length(Vtime))
    disp(['Maximum iterations reached:',num2str(Aero.Wind.V0(3))]);
end
    
figure(1312)
clf
plot(Vtime(1:idT),Power(1:idT))

