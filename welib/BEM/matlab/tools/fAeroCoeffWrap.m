function [Cl Cd Cn Ct CnForAI CtForTI fs] = fAeroCoeffWrap(e,alpha,phi,chord,Vrel,Re,Fperf,WT,Algo)
fs=0;
if(e==0)
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
   e=1:length(Vrel);  %<-------- this will crash, you need the nrom of Vrel if Vrel starts to be the unsteady one
else
    Vreln=norm(Vrel);
end



if(Algo.bPrescribedGamma)
    Gamma=interp1(Algo.r_PrescribedGamma,Algo.PrescribedGamma,WT.Rotor.r(e),'pchip','extrap');
    Cl=2*Gamma./(Vreln.*chord(e));
    Cd=Cl/(Algo.ClOverCd);
else
    bThicknessInterp=logical(1);
% keyboard
    if(Algo.bCl2piAlpha)
        Cl=2*pi*sind(alpha);
        Cd=Cl/(Algo.ClOverCd);%/((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))^(1/20);
    elseif(Algo.bCl2piAlphaCorr)
        Cl=2*pi*sind(alpha).*(1-exp(-(WT.Rotor.r(e)-WT.Rotor.rhub)/(WT.Rotor.R-WT.Rotor.rhub)*1/0.1 ));
        %;*acos(exp(-((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))));
        Cd=Cl/(Algo.ClOverCd);%/((Rotor.r(e)-Rotor.rhub)/(Rotor.R-Rotor.rhub))^(1/20);
    else
        if(~isequal(WT.Sources.Format,'flex')) % bHawc, Hawc2 format with Pc files
            ClCdCm= fAeroCoeff(alpha,WT.Profiles,WT.Rotor.ProfileSet(:,e),WT.Rotor.thickness_rel(e),Re,Algo.bReInterp,bThicknessInterp,Algo.bRoughProfiles);
            Cl=ClCdCm(:,1)';
            Cd=ClCdCm(:,2)';
        else
            ne=length(e);
            Cl=zeros(1,ne);
            Cd=zeros(1,ne);
            for ie=1:ne
                % to be done
                ee=WT.Rotor.ProfileSet(2,e(ie));
                % Badly programmed, what if all the alphas are not the same,
                % then the use of a table is bad
                try
                    Cd(ie)= interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cd(:,ee)  , alpha(ie));
                catch
                    kbd
                end
                if(Algo.bDynaStall)
                    % dynamic stall
                    % interpolation from data
                    f_st=interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.f_st(:,ee) , alpha(ie));
                    Clinv=interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cl_inv(:,ee) , alpha(ie));
                    Clfs= interp1(WT.Profiles.alpha(:,ee) , WT.Profiles.Cl_fs(:,ee)  , alpha(ie));
                    % dynamic stall model
                    tau=4 * WT.Rotor.chord(e) / Vreln;            
                    fs=f_st + ( WT.Aero.last.fs-f_st )*exp(- Algo.dt / tau);
                    Cl(ie)=fs*Clinv+(1-f_st)*Clfs;   
                else
                    % static aerodynamic coefficients
                    Cl(ie)= interp1(WT.Profiles.alpha(:,ee) ,WT.Profiles.Cl(:,ee)  , alpha(ie));
                end    
            end
        end
    end
end
if(Algo.bNoDrag)
    Cd=0*Cd;
end


% performance correction on Cl
if(isfield(Algo.BEM,'bTipLossCl') && Algo.BEM.bTipLossCl)
    if(isequal(Algo.BEM.TipLossMethod,'Shen'))
        warning('Shen and TipLoss Cl are inconsistent... heu, Im not sure...')
    else
        FperfCl=2/pi*acos(exp(-63*(1-WT.Rotor.r(e)/WT.Rotor.R)));
        Cl=Cl.*FperfCl;
    end
end

% Normal and tangential
CnNoDrag=Cl.*cosd(phi);
CtNoDrag=Cl.*sind(phi);
% THIS IS FOR WIND TURBINE, opposite for propellers..
Cn=Cl.*cosd(phi)+Cd.*sind(phi);
Ct=Cl.*sind(phi)-Cd.*cosd(phi);
% performance correction on Cn Ct
Cn=Fperf.*Cn;
Ct=Fperf.*Ct;
if(Algo.bAIDrag)
    CnForAI=Cn;
else
    CnForAI=Fperf.*CnNoDrag;
end
if(Algo.bTIDrag)
    CtForTI=Ct;
else
    CtForTI=CtNoDrag;
end



