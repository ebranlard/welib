function [Cd]=fCdFromCf(xc, Cf, I_sup,I_inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cl from Cp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(xc);
x_up =xc(I_sup);
Cf_up=Cf(I_sup);
x_inf =xc(I_inf);
Cf_inf=Cf(I_inf);
Cd=trapz(x_up,Cf_up)+trapz(x_inf,Cf_inf); % Note, plus sign because all indexes are ascending
%Cd=polyarea(xc,Cf)
end