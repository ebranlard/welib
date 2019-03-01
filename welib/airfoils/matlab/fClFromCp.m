function [Cl]=fClFromCp(xc, Cp, I_sup,I_inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cl from Cp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(xc);
x_up =xc(I_sup);
Cp_up=Cp(I_sup);
x_inf =xc(I_inf);
Cp_inf=Cp(I_inf);
Cl=trapz(x_up,Cp_up)+trapz(x_inf,Cp_inf); % Note, plus sign because all indexes are ascending
%Cl=polyarea(xc,Cp);

end