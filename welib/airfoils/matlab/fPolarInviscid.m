function [Cl_inv_2pi Cl_inv_2pi_sin alpha0]=fPolarInviscid(alpha,Cl)
if max(abs(alpha))>3.2
    % Clearly we are in degrees
    alpha_rad=alpha*pi/180;
else
    % We are most likely in radians
    alpha_rad=alpha;
end

[alpha0, ~, alpha0_rad]=fAlpha0(alpha,Cl);

slope=fAlphaSlope(alpha0_rad,alpha_rad,Cl);

% if abs(slope-2*pi)>pi/4
%     keyboard
%     error('Slope far from 2pi!')
% end
Cl_inv_2pi     = slope*(alpha_rad-alpha0_rad);
Cl_inv_2pi_sin = slope*sin(alpha_rad-alpha0_rad);



