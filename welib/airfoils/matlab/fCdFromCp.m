function [Cd]=fCdFromCp(yc, Cp, I_sup,I_inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cl from Cp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(yc);
y_up =yc(I_sup);
Cp_up=Cp(I_sup);
y_inf =yc(I_inf);
Cp_inf=Cp(I_inf);
Cd=trapz(y_up,Cp_up)+trapz(y_inf,Cp_inf); % Note, plus sign because all indexes are ascending
%Cd=polyarea(yc,Cp);
end