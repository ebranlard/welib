function [R WT] = fBEMunsteady(WT,Sim,Wind,Algo)
%% Initializations
nt=length(Algo.Vtime);
Omega=Sim.Run.Omega;% Used to use WT.Rotor.Omega 
nB   =WT.Rotor.nB;
ne   =WT.Rotor.ne;
dt=Algo.dt;

psi=0; %[rad] !!!

% Forces and Power 
BladeTorque=zeros(nB,nt);
BladeThurst=zeros(nB,nt);
Flap=zeros(nB,nt);
Edge=zeros(nB,nt);

Torque=zeros(1,nt);
Thrust=zeros(1,nt);
Power =zeros(1,nt);
Pitch =zeros(1,nt);
HubV  =zeros(1,nt);


% Results
R.CP=0;
R.CT=0;
R.a= 0;
R.aprime=0;
R.Pn=0;
R.Pt=0;
R.Cl=0;
R.Cd=0;
R.Gamma=0;
R.alpha=0;
R.CTloc=0;
R.CQloc=0;

if(Algo.bUnsteadyStorage)
    R.psi=zeros(nt,nB);
    R.Pn =zeros(nt,nB,ne);
    R.Pt =zeros(nt,nB,ne);
    R.Cl =zeros(nt,nB,ne);
    R.Cd =zeros(nt,nB,ne);
    R.Gamma =zeros(nt,nB,ne);
    R.alpha =zeros(nt,nB,ne);
    R.Un =zeros(nt,nB,ne);
    R.Ut =zeros(nt,nB,ne);
    Vpsi0=mod(0:(360/nB):(360/nB)*(nB-1),360);
%     MKhi=zeros(nt,nB);
%     MW  =zeros(nt,3,nB,ne);
%     MW0 =zeros(nt,3,nB,ne);
end



%% Time loop
for idT=1:nt
    t=Algo.Vtime(idT);
%     if(abs(mod(t,2*pi/Omega))<=dt/2)
    if(abs(mod(psi,2*pi))<Omega*dt)
        disp(['Turn ' num2str(round(t/(2*pi/Omega))) ' - t=' num2str(t)])
    end
    
    %%% Time dependent parameters
    WT.Controller.pitch=WT.Controller.fpitch(t);
    WT.Nacelle.yaw=WT.Controller.fyaw(t);
    Wind.V0=Algo.VV0_of_t(:,mod(idT,length(Algo.VV0_of_t(1,:)))+1);

    x=[0 psi];    %psi [rad]
    v=[0 Omega];  %Omega [rad/s]
    [BEM WT]=fBEM(x,v,1, WT,Sim,Wind,Algo);

    
       
   
    %%%% Time Parameters storage
    Pitch(idT) = WT.Controller.pitch;
    Torque(idT) = WT.Aero.Torque;
    Thrust(idT) = WT.Aero.Thrust;
    Power(idT)  = WT.Aero.Power;
    HubV(idT)=Wind.V0(3);
    %MG(idT)  = Generator.fMoment(v(2));    
    Flap(:,idT)  = WT.Aero.Flap;       
    Edge(:,idT)  = WT.Aero.Edge;       
    
    if(Algo.bUnsteadyStorage)
        R.Cl(idT,1:nB,1:ne)    = BEM.Cl    ; 
        R.Cd(idT,1:nB,1:ne)    = BEM.Cd    ; 
        R.Pn(idT,1:nB,1:ne)    = BEM.Pn    ; 
        R.Pt(idT,1:nB,1:ne)    = BEM.Pt    ; 
        R.Un(idT,1:nB,1:ne)    = BEM.Un    ; 
        R.Ut(idT,1:nB,1:ne)    = BEM.Ut    ; 
        R.Gamma(idT,1:nB,1:ne) = BEM.Gamma ; 
        R.alpha(idT,1:nB,1:ne) = BEM.alpha ; 
        R.psi(idT,1:nB) = mod(Vpsi0 + psi*180/pi,360); % [deg]
    end
    
    %%% Rotation of the blades
    psi=mod(psi+Omega*dt,2*pi) ; %[rad]

end % time loop
R.Power2=WT.Rotor.Omega*Torque;

R.Pitch={Pitch};
R.Power={Power};
R.Torque={Torque};
R.Thrust={Thrust};
R.Flap={Flap};
R.Edge={Edge};



