function qPrime=fqPrime(t,q)

include_globals;

x=q(1); theta=q(2); xDot=q(3); thetaDot=q(4);

if bNoForces
    FHydro=0;
    FWind=0;
    TauHydro=0;
    TauWind=0;
    eta=0;
else
                                            %     [eta, FHydro,TauHydro, u , dudt,vk,fdrag,finertia,ftot] = fHydroCalc(t,q,vf,h,zBot,g,DSpar,vA,vphases,nz,Cm,CD,rhow);
    [eta, FHydro, TauHydro, ~,~, ~, Stored_u, Stored_dudt, dFtot, dMtot, dFdrag, dFinertia]=  fHydroCalcFinal(t,q,vf,vk,vA,vphases,g,rhow,h,zBot,DSpar,Cm,CD,nz,bWheeler,bFloating,bVzTo0);
    if bNoWind
        FWind=0;
        TauWind=0;
    else
        Uw=interp1(vt,vUw,t);
        [V, FWind,TauWind,Ct]  = fWindCalc(t,q,Uw,zHub)  ; 
    end
end
rhs=M\(-B*[xDot; thetaDot] -C*[x; theta] + [FHydro + FWind; TauHydro + TauWind] );
qPrime=[xDot; thetaDot; rhs];


%% Additional variables
if length(q)>4
    qPrime=[qPrime; -gamma_CT*(q(5)-Ct)];
end
%% Some storage of time dependent variables
% it=floor(t/dt)+1;
% vFHydro(end+1)   = FHydro   ; 
% vFWind(end+1)    = FWind   ; 
% vTauHydro(end+1) = TauHydro ; 
% vTauWind(end+1)  = TauWind ; 
% veta(end+1)=eta;
% vtbis(end+1)=t;
return

