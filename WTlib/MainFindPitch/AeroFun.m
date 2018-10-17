
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AeroFun functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,par] = AeroFun(beta,par,Data,flag)

res = 1;
while res > par.res;
    % finding wind speed and angle of attack
    Ut = par.Omega*par.s.*(1+par.at);
    Un = par.Uinf*(1-par.an);
    U = sqrt(Un.^2+Ut.^2);
    psi = atan2(Un,Ut); %(this returns at positive flow angle)

    alpha = psi+par.theta-beta;
    ClCdCm = LiftDataFun(alpha*180/pi,par,Data);
    Ct = ClCdCm(:,1).*sin(psi)-ClCdCm(:,2).*cos(psi);
    Cn = ClCdCm(:,1).*cos(psi)+ClCdCm(:,2).*sin(psi);

    % Induced velocities
    an_old = par.an;
    at_old = par.at;
    if par.TipLoss == 1 % Prandtl's tip loss, as described in M.O.L.Hansen's book
        Ftiploss = Data.Nb*(par.sfull(end)-par.s)./(2*par.s.*sin(psi));
        Ftiploss = 2./pi*acos(exp(-Ftiploss));
    else
        Ftiploss = 1;
    end
    if par.Induction == 1 % induction as done in HAWC2
        CT = (U.^2.*Cn.*par.c*Data.Nb)./(2*pi.*par.s*par.Uinf^2);
        CT = CT./Ftiploss;
        temp = par.k(4)*CT.^3+par.k(3)*CT.^2+par.k(2)*CT+par.k(1);
        par.an = par.relax*par.an+(1-par.relax)*temp;
        par.at = (U.^2.*Ct.*par.c*Data.Nb)./(8*pi.*par.s.^2.*(1-par.an)*par.Uinf*par.Omega);
    end
    res = norm(an_old-par.an)+norm(at_old-par.at);
end
Ft = 1/2*par.rho*par.c.*U.^2.*Ct;
Fn = 1/2*par.rho*par.c.*U.^2.*Cn;

if par.r > 0
    P = trapz(par.sfull',[Data.Nb*Ft*par.Omega;0].*(par.sfull'));
else % skip center point if no hub element
    P = trapz(par.sfull',[0;Data.Nb*Ft*par.Omega;0].*(par.sfull'));
end
if nargin == 4
    P = P-par.Pref;
end
end

