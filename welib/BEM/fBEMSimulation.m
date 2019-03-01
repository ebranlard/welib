function [ BEMout  ] = fBEMSimulation( input_args )


tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Algo.Steady)
    nSimu=size(Simulation.CombinedCase.Cases,1);
    if(nSimu>0)
        disp(['Performing Combined Case'])
        Cases=Simulation.CombinedCase.Cases;
        Thrust=zeros(1,nSimu);
        Flap=zeros(1,nSimu);
        Edge=zeros(1,nSimu);
        Power=zeros(1,nSimu);
        CP=zeros(1,nSimu);
        CT=zeros(1,nSimu);
        vA=zeros(nSimu,length(Rotor.r));
        vAprime=zeros(nSimu,length(Rotor.r));
        for iSim=1:nSimu
            Aero.Wind.V0=[0;0;Cases(iSim,1)];
            Rotor.Omega=Cases(iSim,2);
            Controller.pitch=Cases(iSim,3);          
            %%%
            if(Algo.Steady && Algo.NumSect>1)
                BEM = fBEMpseudo_steady();
            else
                BEM = fBEMsteady();
            end
            %%%
            vBEM(iSim)=BEM;
            BEMout.Thrust(iSim)= BEM.Thrust;
            BEMout.Power(iSim)= BEM.Power;
            BEMout.Edge(iSim)= BEM.Edge;
            BEMout.Flap(iSim)= BEM.Flap;
            BEMout.CP(iSim)= BEM.CP;
            BEMout.CT(iSim)= BEM.CT;
            BEMout.vA(iSim,:)= BEM.a;
            BEMout.vAprime(iSim,:)= BEM.aprime;
        end
    else
        disp(['Performing Parametric Study'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Parametric Study
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Param=Simulation.Parametric;
        if(Param.Pitch(3)==0)
            Vpitch=Param.Pitch(1);
        else
            Vpitch=Param.Pitch(1):Param.Pitch(3):Param.Pitch(2);
        end
        if(Param.Omega(3)==0)
            Vomega=Param.Omega(1);
        else
            Vomega=Param.Omega(1):Param.Omega(3):Param.Omega(2);
        end
        if(Param.WS(3)==0)
            Vws=Param.WS(1);
        else
            Vws=Param.WS(1):Param.WS(3):Param.WS(2);
        end
            
                        
        Thrust=zeros(length(Vpitch),length(Vomega),length(Vws));
        Power=zeros(length(Vpitch),length(Vomega),length(Vws));
        Flap=zeros(length(Vpitch),length(Vomega),length(Vws));
        Edge=zeros(length(Vpitch),length(Vomega),length(Vws));       
        for ip=1:length(Vpitch)
            Controller.pitch=Vpitch(ip); 
            disp(['Pitch:',num2str(Vpitch(ip))]);
            for io=1:length(Vomega)
                Rotor.Omega=Vomega(io);
                disp(['   Omega:',num2str(Vomega(io))]);
                for iws=1:length(Vws)
                    Aero.Wind.V0=[0;0;Vws(iws)];
                    %%%
                    if(Algo.Steady && Algo.NumSect>1)
                        BEM = fBEMpseudo_steady();
                    else
                        BEM = fBEMsteady();
                    end
                    %%%
                    Thrust(ip,io,iws)= BEM.Thrust;
                    Power(ip,io,iws)= BEM.Power;
                    Edge(ip,io,iws)= BEM.Edge;
                    Flap(ip,io,iws)= BEM.Flap;
                end
            end
        end
    end
    
    
else
    if(Algo.DOF==0)
        UnsteadyBEM
    else
        Runge
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc()

if(Simulation.KiloOut)
    Power=Power/1000;
    Thrust=Thrust/1000;
    Flap=Flap/1000;
    Edge=Edge/1000;
  
end

end

