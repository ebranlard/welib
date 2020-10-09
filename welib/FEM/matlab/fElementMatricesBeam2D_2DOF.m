function [ke,me]=fElementMatricesBeam2D_2DOF(EI,L,Mass,MMFormulation)
% Stiffness and mass matrices for Hermitian beam element with 2DOF per node (u theta)
%
% INPUTS
%     EI=E*I: E=elastic modulus, I=second moment of inertia of cross-section
%     L :    Element length
%     M :    Element mass = rho * A * L [kg]
%     MMFormulation = 1 - consistent mass matrix
%                   = 2 - lumped mass matrix
%                   = 3 - diagonal mass matrix
% OUTPUTS
%     ke: Element stiffness matrix (4x4)   
%     me: Element mass matrix (4x4)
%
% AUTHOR: E.Branlard
%
% --- Stiffness matrix
ke=EI/(L^3)*[12      6*L   -12      6*L ; ...
             6*L   4*L^2   -6*L   2*L^2 ; ...
            -12     -6*L    12     -6*L ; ...
             6*L   2*L^2   -6*L   4*L^2];
% ---  Mass matrix
if MMFormulation==1
    % Consistent Formulation
    me=Mass/420*[156    22*L   54    -13*L  ; ...
                 22*L  4*L^2  13*L  -3*L^2  ; ...
                 54    13*L   156      -22*L; ...
                -13*L -3*L^2 -22*L   4*L^2] ; 
elseif MMFormulation==2
    % Lumped formulation
    me=diag([Mass/2  0  Mass/2  0]); % TODO TODO MAY BE WRONG
else
    % Diagonal formulation
    me=Mass*diag([1/2  L^2/78  1/2  L^2/78]);
%     alpha = 17.5;
%     mBeam = rho*A*L / 2 * ...  % lumped
%         [ 1   0              0  0
%     0   alpha*L^2/210  0  0
%     0   0              1  0
%     0   0              0  alpha*L^2/210 ];
%     mElement = mBeam;


end


