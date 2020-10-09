function [ke,me]=fElementMatricesFrame3DOF(E,I,L,A,rho,theta,MMFormulation)
% Stiffness and mass matrices for a 2D frame element: (u1 v1 theta1 u2 v2 theta2)
% Elongation (u)  and bending (v,theta)
%
%  INUPUTS
%     k - element stiffness matrix (size of 6x6)   
%     m - element mass matrix (size of 6x6)
%     E - elastic modulus 
%     I - second moment of inertia of cross-section
%     L - element length
%     A - Area of beam cross-section
%     rho - mass density (mass per unit volume)
%     theta - angle between the local and global axes                           
%            is positive if the local axis is in the ccw direction from
%            the global axis
%     MMFormulation = 1 - consistent mass matrix
%                   = 2 - lumped mass matrix
%                   = 3 - diagonal mass matrix
%
% AUTHOR: E.Branlard, inspired by frame2d
%

% stiffness matrix at the local axis
a=E*A/L;
c=E*I/(L^3);
kl=[
    a   0      0        -a    0       0      ; ...
    0   12*c   6*L*c    0   -12*c    6*L*c   ; ...
    0   6*L*c  4*L^2*c  0   -6*L*c  2*L^2*c  ; ...
    -a  0      0        a     0       0      ; ...
    0   -12*c  -6*L*c   0   12*c    -6*L*c   ; ...
    0   6*L*c  2*L^2*c  0   -6*L*c  4*L^2*c] ; 
% rotation matrix
r=[ cos(theta)  sin(theta)  0   0           0           0;...
   -sin(theta)  cos(theta)  0   0           0           0;...
    0          0            1   0           0           0;...
    0          0            0   cos(theta)  sin(theta)  0;...
    0          0            0  -sin(theta)  cos(theta)  0;...
    0          0            0   0           0           1];

% --- Mass Matrix
mass=rho*A*L;
if MMFormulation==1
    % consistent mass matrix
    mm=mass/420;
    ma=mass/6;
    ml=[2*ma  0           0         ma      0              0    ;...
        0     156*mm      22*L*mm   0       54*mm       -13*L*mm;...
        0     22*L*mm     4*L^2*mm  0       13*L*mm    -3*L^2*mm;...
        ma    0           0        2*ma       0            0    ;...
        0     54*mm       13*L*mm   0       156*mm      -22*L*mm;...
        0    -13*L*mm    -3*L^2*mm  0      -22*L*mm    4*L^2*mm];
elseif MMFormulation==2
    % lumped mass matrix
    ml=mass*diag([0.5  0.5  0  0.5  0.5  0]);
else
    % diagonal mass matrix
    ml=mass*diag([0.5  0.5  L^2/78  0.5  0.5  L^2/78]);
end

% Conversion to global system
me = r'*ml*r;
ke = r'*kl*r;
