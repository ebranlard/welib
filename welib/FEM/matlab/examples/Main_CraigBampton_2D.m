%% Documentation   
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

%% Main Parameters
nel      = 10;
nModesCB = 6 ;
BC       = 'clamped-free';  % Boundary condition: free-free or clamped-free
BC       = 'free-free';  % Boundary condition: free-free or clamped-free

% --- Algo specific params
ndof  = 2 ; % Degrees of Freedom per Node: 2 for 2D beams
% --- Tube Beam - Structural data
L   = 100      ;
EI0 = 1.654e+12; 
m0  = 1.026e+04;
A   = 1        ;
rho = m0       ;

%% --- Building mass and stiffness matrix
gravity = 9.81;
Mtop    = 0   ;
[MM,KK,KKg,x]=fBeamMatrices2D_2DOF(L,EI0,m0,nel,gravity,Mtop,[],[],true,true);
% KK=KK+KKg;
fprintf('Number of DOF of full system:    %d\n',size(KK,1));

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
U_full=Q_full(1:ndof:end,:);
V_full=Q_full(2:ndof:end,:);

% nModesPlot=4;
% for i=1:nModesPlot
%     figure,clf,hold all,box on;
%     fact=1 / U_full(end,i);
%     plot(x,U_full(:,i)*fact)
%     plot(x,V_full(:,i)*fact)
%     legend('U','V')
%     title(sprintf('Full system mode %d', i))
% end


%% --- CB reduction
Im = [-ndof+1:0]+size(MM,1);  % Selecting top node (last node) as master 
xs = x_nobc(1:end-1);

[fc,Mr,Kr,~,~,Psi_sq, ~, ~, Psi_sm] = fCraigBamptonReduce(Im,nModesCB,MM,KK);
Psi_sm(abs(Psi_sm)<1e-8)=0;

fprintf('fc=%8.3f  (first constrained frequency)\n',fc);
% [fc,Mr,Kr] = fGuyanReduce(Im,MM,KK);
fprintf('Number of DOF after reduction:  %d\n',size(Kr,1))

% --- Guyan Modes
U_Guy=Psi_sm(1:ndof:end,:);
V_Guy=Psi_sm(2:ndof:end,:);
figure
for i=1:length(Im)
    subplot(1,length(Im),i); hold all; box on
    plot(xs,U_Guy(:,i))
    plot(xs,V_Guy(:,i))
    legend('U','V')
    title(sprintf('Guyan mode %d',i))
end
% --- CB modes
U_CB=Psi_sq(1:ndof:end,:);
V_CB=Psi_sq(2:ndof:end,:);
figure
for i=1:floor(nModesCB/2)
    subplot(1,floor(nModesCB/2),i); hold all; box on
    plot(xs,U_CB(:,i))
    plot(xs,V_CB(:,i))
    legend('U','V')
    title(sprintf('CB mode %d',i))
end

%% --- EVA of reduced system
[Q,Lambda]=eig(Kr,Mr);
Omega2=diag(Lambda);
[Omega2,Isort]=sort(Omega2);
Q=Q(:,Isort);
f= sqrt(Omega2)/(2*pi);
for i=1:min(size(Kr,1),8);
    fprintf('                                             %.3f \n',f(i));
end

