function [M K C]=getMatrices()
global Tower Shaft Rotor Nacelle Generator Algo

M_tot=Rotor.M+Nacelle.M;

% stiffness matrix
K=[Tower.k 0 0; 0 0 0 ;  0 0 Shaft.k];
% Additional damping
C=[0 0 0
    0 0 0
    0 0 Algo.damp];

% maxx matrix
M=[M_tot  0 0 ; ...
         0 Rotor.I + Generator.I  Generator.I  ;...
         0  Generator.I Generator.I];
end