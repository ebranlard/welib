% function [ a aprime CT] = fInductionCoefficients(a_last,F,Ftip,CnForAI,Ct,lambda_r,sigma,phi,Algo)
function [ a aprime CT_loc] = fInductionCoefficients(a_last,Vrel,Un,Ut,V0_in3,V0_in4,nW_in4,omega,chord,F,Ftip,CnForAI,CtForTI,lambda_r,sigma,phi,Algo) ;
if(length(lambda_r)>1)
    n=size(Vrel,2);
    if(size(Vrel,1)==3)
       %Vrel is a vector of three component velocity
       Vreln=zeros(1,n);
       for i=1:n
           Vreln(i)=norm(Vrel(:,i));
       end
   else
      Vreln=Vrel;
  end 
else
    Vreln=norm(Vrel);
end





%normal expression, the default one
a=1./( (4*F.*sind(phi).^2) ./(sigma.*CnForAI)+1 );
% Thrust coefficient from the momentum theory => alast
% CT=(1-a_last).^2.*sigma.*CnForAI./((sind(phi)).^2)
CT_loc=Vreln.^2.*sigma.*CnForAI./norm(V0_in3)^2    ; % that's a CT loc
% keyboard
%hawc CT CT = (U.^2.*sigma)./(par.Uinf^2);
%         CT = (U.^2.*Cn.*par.c*Data.Nb)./(2*pi.*par.s*par.Uinf^2);
%         CT = CT./Ftiploss;

if(isfield(Algo.BEM,'bTipLossCl') && Algo.BEM.bTipLossCl)% if(Algo.bTipLoss && isequal(Algo.TipLossMethod,'Shen')) 
    Y1=4*Ftip.*sind(phi).^2./(sigma.*CnForAI);
    Y2=4*Ftip.*sind(phi).*cosd(phi)./(sigma.*CtForTI);
    a=(2+Y1-sqrt(4*Y1.*(1-Ftip)+Y1.^2) ) ./(2*(1+Ftip.*Y1));
end

%%% Correction for high induction, BEM breakdown 
if(isequal(Algo.BEM.Kind,'Bramwell')) % will be improved when I'll know more about yaw, tilt models
    ac=0.2;
    if a_last<=ac
        fg=1;
    elseif a_last>ac
        fg=ac/a_last*(2-ac/a_last);
    end
%     W_z_qs=-nB*BEM.L*cosd(phi)/(4*pi*rho*r(e)*BEM.F*norm(V0_in4+fg*nnW_in4 ));
%     W_y_qs=-nB*BEM.L*sind(phi)/(4*pi*rho*r(e)*BEM.F*norm(V0_in4+fg*nnW_in4 ));
    W_z_qs=Vreln.^2.*CnForAI/(4*F*norm(V0_in4+fg*nW_in4 ) )*sigma;
    W_y_qs=Vreln.^2.*CtForTI/(4*F*norm(V0_in4+fg*nW_in4 ) )*sigma;
    a=W_z_qs/norm(V0_in3);
    aprime=W_y_qs/(lambda_r*norm(V0_in3));
    if(Algo.bSteady)
        %relaxation
        a=a*Algo.relaxation+(1-Algo.relaxation)*a_last;
    end
else
    
    % --------------------------------------------------------------------------------
    % --- Hight thrust correction 
    % --------------------------------------------------------------------------------
    [ a CT_loc] = fCorrectionHighThrust(Algo.BEM.CTCorrection,a, CnForAI ,phi, a_last,sigma, F,Ftip, CT_loc );

    if(Algo.bSteady)
        %relaxation
        a=a*Algo.relaxation+(1-Algo.relaxation)*a_last;
    end

    %%% Swirl
    if(Algo.bSwirl)
        if(Algo.BEM.bTipLoss && isequal(Algo.BEM.TipLossMethod,'Shen') && isfield(Algo.BEM,'bTipLossCl') && Algo.BEM.bTipLossCl)
            aprime=1./((1-a.*Ftip).*Y2./ ((1-a)-1));  % !!!!!!!!!!!!!!!!!!!!!! doubts on this one
        else
            aprime=1./( (4*F.*sind(phi).*cosd(phi)) ./(sigma.*CtForTI)  -1 );
            if(isequal(Algo.BEM.SwirlMethod,'AeroDyn'))
                SwlAng=1+4*a.*F.*(1-a)./lambda_r.^2;
                aprime=0.5*(sqrt(SwlAng)-1);
            elseif(isequal(Algo.BEM.SwirlMethod,'Hawc'))
                aprime=(Vreln.^2.*CtForTI.*sigma)./(4*(1-a)*norm(V0_in4).^2.*lambda_r);
                %            aprime=(Vrel^2*Ct*chord*B)/(8*pi*Algo.r^2*(1-a)*V0*Algo.Omega);
                %             fprintf('%.7f \n',abs(aprime0-aprime));
                %         par.at = (U.^2.*Ct.*par.c*Data.Nb)./(8*pi.*par.s.^2.*(1-par.an)*par.Uinf*par.Omega);
            end

        end
    else
        aprime=a*0;
    end
    if(sum(isnan(a))>0 )
        warning('BEM is crashing')
        keyboard
        error('Break')
    end


end
