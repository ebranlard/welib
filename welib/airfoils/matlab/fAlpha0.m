function [alpha0 alpha0_deg alpha0_rad]=fAlpha0(alpha,Cl)
% returns the angle of attach of zero lift in the same dimension as the input (degree or radian)
if max(abs(alpha))>3.2
    % Clearly we are in degrees

    % find alpha_0 between -20 and +20 (can be improved)
    Im20=whichvalue(alpha,-20);
    Ip20=whichvalue(alpha,20);
    I=Im20:Ip20;
    ialpha0=I(find(Cl(I)>0,1,'first')-1);
    alpha0=interp1(Cl(ialpha0:ialpha0+1),alpha(ialpha0:ialpha0+1),0);
    alpha0_deg=alpha0;
    alpha0_rad=alpha0*pi/180;
else 
    % We are most likely in radians

    % find alpha_0 between -20 and +20 (can be improved)
    Im20=whichvalue(alpha,-20*pi/180);
    Ip20=whichvalue(alpha,20*pi/180);
    I=Im20:Ip20;
    ialpha0=I(find(Cl(I)>0,1,'first')-1);
    alpha0=interp1(Cl(ialpha0:ialpha0+1),alpha(ialpha0:ialpha0+1),0);
    alpha0_rad=alpha0;
    alpha0_deg=alpha0*180/pi;
end
