function [Cl_inv_2pi Cl_inv_2pi_sin alpha0]=fPolarInviscid(alpha,Cl)
if max(abs(alpha))>3.2
    % Clearly we are in degrees
    alpha_rad=alpha*pi/180;
else
    % We are most likely in radians
    alpha_rad=alpha;
end

[alpha0, ~, alpha0_rad]=fAlpha0(alpha,Cl);

Cl_inv_2pi     = 2*pi*(alpha_rad-alpha0_rad);
Cl_inv_2pi_sin = 2*pi*sin(alpha_rad-alpha0_rad);



