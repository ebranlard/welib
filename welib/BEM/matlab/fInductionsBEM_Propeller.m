function [ uin uit CT_loc] = fInductionsBEM_Propeller(uin_prev,uit_prev,Un,Ut,V0,omega,chord,F,Ftip,CnForAI,CtForTI,r,sigma,phi,Cl,Algo) ;

% --------------------------------------------------------------------------------
% --- Case with no free stream 
% --------------------------------------------------------------------------------
if norm(V0)< 1e-6
    if isequal(Algo.BEM.Kind,'FormulationCl')
        % --------------------------------------------------------------------------------
        % --- Formulation using Cl
        % --------------------------------------------------------------------------------
        if(Algo.bSwirl)
            K  = sigma./(4*F).*Cl;

            uit = omega*r .*K.^2 /2* ( sqrt(1+4./K.^2)-1);
            % --- relaxation
            uit = uit*Algo.relaxation+(1-Algo.relaxation)*uit_prev;

            uin = sqrt(uit.*(omega*r-uit));
        else
            error('Propeller without swirl TODO')

        end
    elseif isequal(Algo.BEM.Kind,'FormulationCnCt')
        % --------------------------------------------------------------------------------
        % --- Formulation using Cn and Ct 
        % --------------------------------------------------------------------------------
        if(Algo.bSwirl)
            %  Buggy:
            Kt  = sigma./(4*F.*sind(phi).*cosd(phi)).*CtForTI;
            uit = omega*r./(1./Kt-1) ;
            %
            uin = CnForAI./CtForTI.*uit;
            % --- relaxation
            uin=uin*Algo.relaxation+(1-Algo.relaxation)*uin_prev;
        else
            Kn  = sigma./(4*F.*sind(phi).*cosd(phi)).*CnForAI;
            uin = Kn.* omega*r;
            uit=0;
            Ibad= abs(phi)<1e-9;
            if sum(Ibad)>0
                fprintf('Position where the flow is zero %d \n',sum(Ibad));
                phi(Ibad)=1e-3;
                Kn  = sigma./(4*F.*sind(phi).*cosd(phi)).*CnForAI;
                uin = Kn.* omega*r;
            end

            % --- relaxation
            uin=uin*Algo.relaxation+(1-Algo.relaxation)*uin_prev;
        end
    else
        error(['Unknwon Algo.BEM.Kind for propeller induction' Algo.BEM.Kind])
    end
   
else
    error('Propeller with free stream todo')

end

CT_loc=NaN;
