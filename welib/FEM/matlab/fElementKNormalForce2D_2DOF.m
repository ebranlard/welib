function Kg=fElementKNormalForce2D_2DOF(L,N)
% Geometrical stiffness matrix based on a normal force (self-weight, top-mass, centrifugal force)
%
% INPUTS
%   L : length of element
%   N : Normal force (negative in compression)
%
% AUTHOR: E. Branlard
%
Kg = N/(30*L) * [ 36    3*L   -36   3*L
                  3*L   4*L^2 -3*L  -L^2
                  -36   -3*L    36  -3*L
                  3*L   -L^2 -3*L  4*L^2 ];
