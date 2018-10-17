function [R]= fWTSimulation(Code,WT,Sim,Wind,Algo)
R.Code=Code;


fprintf('Code : %s - ',Code)
tic()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCases=size(Sim.CombinedCases,1);
if(nCases>0)
    disp(['Performing Combined Case  '])
    Cases=Sim.CombinedCases;
    R.WS=Cases(:,1);
    R.RPM=Cases(:,2);
    R.PITCH=Cases(:,3);
    R.YAW=Cases(:,4);

    R.Omega=Cases(:,2)*2*pi/60;
    R.lambda=R.Omega*WT.Rotor.R./R.WS;;

    if(Algo.bSteady)
      R.Thrust=zeros(1,nCases);
      R.Flap=zeros(1,nCases);
      R.Edge=zeros(1,nCases);
      R.Power=zeros(1,nCases);
      R.CP=zeros(1,nCases);
      R.CT=zeros(1,nCases);

      R.vA=zeros(nCases,length(WT.Rotor.r));
      R.vAprime=zeros(nCases,length(WT.Rotor.r));

      R.vPn=zeros(nCases,length(WT.Rotor.r));
      R.vPt=zeros(nCases,length(WT.Rotor.r));
      R.vCl=zeros(nCases,length(WT.Rotor.r));
      R.vCd=zeros(nCases,length(WT.Rotor.r));
      R.vAlpha=zeros(nCases,length(WT.Rotor.r));
      R.vCTloc=zeros(nCases,length(WT.Rotor.r));
      R.vCQloc=zeros(nCases,length(WT.Rotor.r));
    end

    for iSim=1:nCases
        fprintf('.');
        Sim.Run=fSetRun(Cases(iSim,:));
        Sim.Run.Omega=Sim.Run.RPM*2*pi/60;
        Sim.Run.lambda=Sim.Run.Omega*WT.Rotor.R/Sim.Run.WS;
        Wind=fSetWind(Wind,Sim);
        Algo.VV0_of_t=[0;0;Sim.Run.WS]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Not pretty
        [ Results ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
        if(Algo.bSteady)
          R.Results(iSim)=Results;
          R.Thrust(iSim)= Results.Thrust;
          R.Power(iSim)= Results.Power;
          R.Edge(iSim)= Results.Edge;
          R.Flap(iSim)= Results.Flap;
          R.CP(iSim)= Results.CP;
          R.CT(iSim)= Results.CT;
          R.vA(iSim,:)= Results.a;
          R.vAprime(iSim,:)= Results.aprime;
          R.vPn(iSim,:)= Results.Pn;
          R.vPt(iSim,:)= Results.Pt;
          R.vCl(iSim,:)= Results.Cl;
          R.vCd(iSim,:)= Results.Cd;
          R.vAlpha(iSim,:)= Results.alpha;
          R.vCTloc(iSim,:)= Results.CTloc;
          R.vCQloc(iSim,:)= Results.CQloc;
        else
          R.Results{iSim} = Results        ; 
          R.Thrust{iSim}  = Results.Thrust ; 
          R.Power{iSim}   = Results.Power  ; 
          R.Edge{iSim}    = Results.Edge   ; 
          R.Flap{iSim}    = Results.Flap   ; 
          R.CP{iSim}      = Results.CP     ; 
          R.CT{iSim}      = Results.CT     ; 
          R.vA{iSim}      = Results.a      ; 
          R.vAprime{iSim} = Results.aprime ; 
          R.vPn{iSim}     = Results.Pn     ; 
          R.vPt{iSim}     = Results.Pt     ; 
          R.vCl{iSim}     = Results.Cl     ; 
          R.vCd{iSim}     = Results.Cd     ; 
          R.vAlpha{iSim}  = Results.alpha  ; 
          R.vCTloc{iSim}  = Results.CTloc  ; 
          R.vCQloc{iSim}  = Results.CQloc  ; 
        end
    end
else
    disp(['Performing Parametric Study'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Parametric Study
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vws=Sim.ParametricWS;
    Vrpm=Sim.ParametricRPM;
    Vpitch=Sim.ParametricPITCH;
    Vyaw=Sim.ParametricYAW;


    R.WS=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));
    R.RPM=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));
    R.PITCH=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));
    R.YAW=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));    
    
    R.lambda=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));    
    R.Omega=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));    

    R.Thrust=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));
    R.Power=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));
    R.Flap=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));
    R.Edge=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));       
    R.CP=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));       
    R.CT=zeros(length(Vyaw),length(Vpitch),length(Vrpm),length(Vws));       
    for iy=1:length(Vyaw)
        fprintf(['Yaw:',num2str(Vyaw(iy))]);
        for ip=1:length(Vpitch)
            fprintf([' Pitch:',num2str(Vpitch(ip))]);
            for io=1:length(Vrpm)
                fprintf([' RPM:',num2str(Vrpm(io)) '  ']);
                for iws=1:length(Vws)
                    fprintf('.');
                    R.WS(iy,ip,io,iws)= Vws(iws);
                    R.RPM(iy,ip,io,iws)= Vrpm(io); 
                    R.PITCH(iy,ip,io,iws)= Vpitch(ip); 
                    R.YAW(iy,ip,io,iws)=  Vyaw(iy);

                    Sim.Run=fSetRun(Vws(iws),Vrpm(io),Vpitch(ip),Vyaw(iy));
                    Sim.Run.Omega=Sim.Run.RPM*2*pi/60;
                    Sim.Run.lambda=Sim.Run.Omega*WT.Rotor.R/Sim.Run.WS;
                    Wind=fSetWind(Wind,Sim);
                    
                    R.Omega(iy,ip,io,iws)=Sim.Run.Omega;
                    R.lambda(iy,ip,io,iws)=Sim.Run.lambda;

                    [ Results ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
                    R.Results(iy,ip,io,iws)=Results;
                    R.Thrust(iy,ip,io,iws)= Results.Thrust;
                    R.Power(iy,ip,io,iws)= Results.Power;
                    R.Edge(iy,ip,io,iws)= Results.Edge;
                    R.Flap(iy,ip,io,iws)= Results.Flap;
                    R.CT(iy,ip,io,iws)= Results.CT;
                    R.CP(iy,ip,io,iws)= Results.CP;


                end
                fprintf('\n')
            end
        end
    end

    R.Results=squeeze(R.Results);

    R.WS=squeeze(R.WS);
    R.RPM=squeeze(R.RPM);
    R.PITCH=squeeze(R.PITCH);
    R.YAW=squeeze(R.YAW);
    
    R.Omega=squeeze(R.Omega);
    R.lambda=squeeze(R.lambda);



    R.Thrust=squeeze(R.Thrust);
    R.Power=squeeze(R.Power);
    R.Edge=squeeze(R.Edge);
    R.Flap=squeeze(R.Flap);
    R.CP=squeeze(R.CP);
    R.CT=squeeze(R.CT);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc()

if(Sim.bKiloOut)
    if(Algo.bSteady)
      R.Power=R.Power/1000;
      R.Thrust=R.Thrust/1000;
      R.Flap=R.Flap/1000;
      R.Edge=R.Edge/1000;
    else
      warning('For now, can''t use bKiloOut with Unsteady')
    end

end
