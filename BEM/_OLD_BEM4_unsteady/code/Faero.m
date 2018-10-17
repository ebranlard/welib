function [F fs]=Faero(a,dt)
% compute aerodynamical force from aerodynamic data
fs=0;
if(isequal(a.model,'dynastall'))
    % dynamic stall
    % interpolation from data
    f_st=interp1(a.Profile.alpha , a.Profile.f_st , a.alpha);
    Clinv=interp1(a.Profile.alpha , a.Profile.Cl_inv , a.alpha);
    Clfs= interp1(a.Profile.alpha , a.Profile.Cl_fs  , a.alpha);
    % dynamic stall model
    tau=4 * a.Profile.cr / a.Vrel;
    fs=f_st + ( a.fs_prev-f_st )*exp(-dt/tau);
    Cl=fs*Clinv+(1-f_st)*Clfs;
else
    % static aerodynamic coefficients
    Cl= interp1(a.Profile.alpha , a.Profile.Cl  , a.alpha);
end

F=0.5* a.rho * a.Vrel^2 * a.Profile.cr * a.Profile.span * Cl * cos(a.phi);

end