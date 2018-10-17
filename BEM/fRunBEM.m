function [ BEM BEM2] = fRunBEM(WT,Sim,Wind,Algo )
BEM2=[];
%% Some verifications due to change in structure of Algo
if isfield(Algo,'bTipLoss') || isfield(Algo,'bTipLoss')
    error('bTipLoss moved to Algo.BEM')
end
if isfield(Algo,'bHubLoss') || isfield(Algo,'bHubLoss')
    error('bHubLoss moved to Algo.BEM')
end
if isfield(Algo,'Correction') || isfield(Algo.BEM,'Correction')
    error('Correction moved to BEM.CTCorrection and BEM.SwirlMethod')
end

%% Some verifications
if(( Wind.nu~=0 || Sim.Run.YAW ~=0 || WT.Nacelle.tilt~=0))
    if(Algo.BEM.NumSect==1)
         disp('!Warning: You should use different sectors for BEM computations')
    end
else
    if(Algo.BEM.NumSect>1)
        disp('!Warning: There is no need for several sectors - Going back to 1 sector')
        Algo.BEM.NumSect=1;
    end
end
if(Algo.BEM.NumSect>1)
    Algo.BEM.NumSect=max(Algo.BEM.NumSect,4);
end

if(norm(Wind.V0)~=Sim.Run.WS)
    error('Wind and Sim wind speeds disagree in norm')
end


%% Unsteady simulations initialization
% Aero
if((Algo.bSteady && Algo.BEM.NumSect>1) || ~Algo.bSteady)
    if(Algo.bSteady)
        disp('Will perform pseudo-unsteady for wake equilibrium')
        Algo.bYawModel=1; %in this if?
    end
    WT.Aero.last.W=ones(3,WT.Rotor.ne,WT.Rotor.nB).*    Algo.BEM.w_guess;      %temporary induced velocity
    WT.Aero.last.W0=ones(3,WT.Rotor.ne,WT.Rotor.nB).*   Algo.BEM.w_guess;     %temporary induced velocity
    WT.Aero.last.W_qs=ones(3,WT.Rotor.ne,WT.Rotor.nB).* Algo.BEM.w_guess;   %temporary quasistatic induced velocity
    WT.Aero.last.W_int=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;  %temporary intermediate induced velocity
    WT.Aero.last.chi=Sim.Run.YAW; % deg
end
%% Running BEM
% tic()fprintf('BEM run...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Algo.bSteady)
    if(Algo.BEM.NumSect>1)
        BEM = fBEMpseudo_steady();
    else
        if(isequal(Algo.BEM.TipLossMethod,'TipLossDB') && Algo.BEM.bTipLoss)
            % first run without tip loss
            %                  Algo.TipLoss=0;
            if(isfield(Algo.BEM,'bTipLossCl') && Algo.BEM.bTipLossCl)
                Algo.BEM.TipLossMethod='Shen';
            else
                Algo.BEM.TipLossMethod='Glauert';
            end
            fprintf('First BEM run...')
            tic()
            BEM = fBEMsteady(WT,Sim,Wind,Algo);    
            toc()
            % then run with tip loss
            fprintf('\nSecond BEM run...')
            Algo.BEM.TipLossMethod='TipLossDB';
            tic()
            disp('remember there is a choice going on here between standalone or not, to be implemented properly');
            %BEM2 = fBEMsteadyTipLoss(WT,Sim,Wind,Algo,BEM);   
            BEM2=fBEMsteadyTipLossStandAlone(WT,Sim,Wind,Algo);
            toc()
        elseif(isequal(Algo.BEM.TipLossMethod,'PrescribedWake') && Algo.BEM.bTipLoss)
            % we can run directly the new BEM
            disp('New PrescribedWake tip-loss')
            BEM=fBEMsteadyTipLossStandAlone(WT,Sim,Wind,Algo);
        else
            % normal BEM
            BEM = fBEMsteady(WT,Sim,Wind,Algo);
        end
    end
else
    if(Algo.DOF==0)
        WT.Aero.last.W=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;      %temporary induced velocity
        WT.Aero.last.W0=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;     %temporary induced velocity
        WT.Aero.last.W_qs=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;   %temporary quasistatic induced velocity
        WT.Aero.last.W_int=ones(3,WT.Rotor.ne,WT.Rotor.nB).*Algo.BEM.w_guess;  %temporary intermediate induced velocity
        if (abs(Sim.Run.YAW)>0) 
            WT.Controller.yaw=Sim.Run.YAW;
            WT.Controller.fyaw=@(t)Sim.Run.YAW;
            fprintf('!Warning: fBEMUnsteady Yaw used to be provided by controller\n');
        end
        if (abs(Sim.Run.PITCH)>0) 
            WT.Controller.fpitch=@(t)Sim.Run.PITCH;
            WT.Controller.pitch=Sim.Run.PITCH;
            fprintf('!Warning: fBEMUnsteady Pitch used to be provided by controller\n');
        end
        BEM=fBEMunsteady(WT,Sim,Wind,Algo);
    else
        Runge
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Done. \t\t\t(%09.5f s) \n',toc())



end

