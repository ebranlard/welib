function [Cl]=fClFromCf(yc, Cf, I_sup,I_inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cl from Cp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(yc);
y_up =yc(I_sup);
Cf_up=Cf(I_sup);
y_inf =yc(I_inf);
Cf_inf=Cf(I_inf);
Cl=trapz(y_up,Cf_up)+trapz(y_inf,Cf_inf); % Note, plus sign because all indexes are ascending
%Cl=polyarea(yc,Cf);
end