%% Documentation   
% NOTE: Frame 3D, main direction is along x
% 
% 1. Setting up a system with free-free condition
%    - Rigid body modes are obtained (2 in 2D), with 0 frequencies
%    - The rigid body modes obtained by the EVA are two modes of constant slopes but opposite sign
%    - A linear combination of those two can be used to create a mode of pure translation and 0 rotation
% 2. Performing CB reduction of the system:
%    - Guyan modes correponds to rigid body translation and rotation of top node
%    - CB modes are constrained at top node, and free on the other side
% 
% See also: 
%      ThinBeam3D
% 
% Contact: E. Branlard 

%% Init
clear all; close all; clc;
restoredefaultpath;
addpath('../');

%% Main Parameters
nel      = 10;
nModesCB = 6 ;
BC       = 'clamped-free';  % Boundary condition: free-free or clamped-free
% BC       = 'free-free';  % Boundary condition: free-free or clamped-free

% --- Algo specific params
ndof  = 6 ; % Degrees of Freedom per Node: 6 for 3D beams
% --- Tube Beam - Structural data
E   = 210e9    ; % Young modulus [Pa] [N/m^2]
G   = 79.3e9   ; % Shear modulus. Steel: 79.3  [Pa] [N/m^2]
L   = 100      ; % Beam Length [m]
EIy0= 1.654e+12; % Planar second moment of area [m^4]
m0  = 1.026e+04; % Mass per length [kg/m]
EIx0= EIy0*2   ; % Polar second moment of area [m^4]
A   = 1.00     ;
rho = m0       ;
Kt  = EIy0/E*10; % Torsion constant [m^4]

%% --- Building mass and stiffness matrix
% Constant beam properties
[MM,KK,x]=fBeamMatrices3D_Frame6DOF(L,EIx0,EIy0,EIy0,E*A,m0,A,E,G,Kt,nel);
fprintf('Number of DOF of full system:   %d\n',size(KK,1));

%% --- Applying BC
if isequal(BC,'clamped-free')
    for i=1:ndof
        KK(1,:)=[]; KK(:,1) = [];
        MM(1,:)=[]; MM(:,1) = [];
    end
    x_nobc=x(2:end);
else
    x_nobc=x;
end
fprintf('Number of DOF after BC:         %d\n',size(KK,1))


%% --- EVA on system with BC
[Q_full,Lambda]=eig(KK,MM);
Omega2=diag(Lambda);
[Omega2,Isort]=sort(Omega2);
Q_full=Q_full(:,Isort);
f= sqrt(Omega2)/(2*pi);
for i=1:min(size(KK),8);
    fprintf('                                             %.3f \n',f(i));
end
% --- Mode shapes
if isequal(BC,'clamped-free')
    % Adding 0s from BC
    Q_full=[zeros(ndof,size(Q_full,2)); Q_full];
end
for j = 1:3
    U_full{j} = Q_full(    j:ndof:end,:);
    V_full{j} = Q_full((3+j):ndof:end,:);
end

% nModesPlot=5;
% for i=1:nModesPlot
%     figure,clf,hold all,box on;
%     for j = 1:3
%         plot(x,U_full{j}(:,i),'-')
%         plot(x,V_full{j}(:,i),'--')
%     end
%     legend('Ux','Vx','Uy','Vy','Uz','Vz')
%     title(sprintf('Full system mode %d', i))
% end
% 
% 
% %% --- CB reduction
% Im = [-ndof+1:0]+size(MM,1);  % Selecting top node (last node) as master 
% xs = x_nobc(1:end-1);
% 
% [fc,Mr,Kr,~,~,Psi_sq, ~, ~, Psi_sm] = fCraigBamptonReduce(Im,nModesCB,MM,KK);
% Psi_sm(abs(Psi_sm)<1e-8)=0;
% 
% fprintf('fc=%8.3f  (first constrained frequency)\n',fc);
% % [fc,Mr,Kr] = fGuyanReduce(Im,MM,KK);
% fprintf('Number of DOF after reduction:  %d\n',size(Kr,1))
% 
% % --- Guyan Modes
% for j = 1:3
%     U_Guy{j}=Psi_sm(  j:ndof:end,:);
%     V_Guy{j}=Psi_sm(3+j:ndof:end,:);
% end
% figure
% for i=1:length(Im)
%     subplot(1,length(Im),i); hold all; box on
%     for j = 1:3
%         plot(xs,U_Guy{j}(:,i),'-')
%         plot(xs,V_Guy{j}(:,i),'--')
%     end
%     legend('Ux','Vx','Uy','Vy','Uz','Vz')
%     title(sprintf('Guyan mode %d',i))
% end
% % --- CB modes
% for j = 1:3
%     U_CB{j}=Psi_sq(  j:ndof:end,:);
%     V_CB{j}=Psi_sq(3+j:ndof:end,:);
% end
% figure
% for i=1:floor(nModesCB/2)
%     subplot(1,floor(nModesCB/2),i); hold all; box on
%     for j = 1:3
%         plot(xs,U_CB{j}(:,i))
%         plot(xs,V_CB{j}(:,i))
%     end
%     legend('Ux','Vx','Uy','Vy','Uz','Vz')
%     title(sprintf('CB mode %d',i))
% end
% 
% %% --- EVA of reduced system
% [Q,Lambda]=eig(Kr,Mr);
% Omega2=diag(Lambda);
% [Omega2,Isort]=sort(Omega2);
% Q=Q(:,Isort);
% f= sqrt(Omega2)/(2*pi);
% for i=1:min(size(Kr,1),8);
%     fprintf('                                             %.3f \n',f(i));
% end
% 
