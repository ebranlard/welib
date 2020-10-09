%% Documentation
% Setup a finite element model of a tower (with constant properties)
% Performs time integration where the tower is excited by a point force at the top.
% Options:
%    - Perform a Craig-Bampton reduction
%    - Use a "soil element", instead of the clamped Boundary Condition. 
%      The soil element 2x2 stiffness matrix is determined using a 2d frame
%        - The soil element can be added to the system matrix
%        - The soil element can be added to the force (more challenging).
%
% Combinations supported:
%  - No CB             - No Soil El.             - clamped-free - Newmark (dt=0.005), RK4 dt=0.000005)
%  -    CB(top)        - No Soil El.             - clamped-free - Newmark (dt=0.005), RK4 dt=0.0005)
%  - No CB             -    Soil El. (Not force) - free-free    - Newmark (dt=0.005), RK4 dt=0.000005)
%  -    CB(top-bottom) -    Soil El. (Not force) - free-free    - Newmark (dt=0.005), RK4 dt=0.000005)
%
% Combinations not yet supported: 
%  - No CB             -    Soil El. (force)     - free-free    - ode23
%   
%
% Contact: E. Branlard 

%% Init
clear all; close all; clc;
restoredefaultpath;
YAMS_PATH =  fullfile(pwd,'../../../');
addpath(genpath([YAMS_PATH 'fem']));
addpath(genpath([YAMS_PATH 'ode']));

global nit; nit=0

%% Main Parameters
% --- Main flags
bPlotModes      = logical(0)    ;
bCBReduction    = logical(0)    ;
bUseSoilElement = logical(0)    ;
BC              = 'clamped-free'; % Boundary condition: free-free or clamped-free
% BC              = 'free-free'   ; % Boundary condition: free-free or clamped-free
% --- Options for CB reduction (Only if bCBReduction)
CBLeader        = 'top'         ; 
% CBLeader        = 'top-bottom'  ;
nModesCB        = 6             ; 
% --- Options for soil (only if bUseSoilElement)
bSoilAsForce    = logical(1)    ; % Insert soil contribution as a force or in the stiffness matrix
dx_soil         = 0.01          ; % If using a "soil" element 
% --- Model 
nel             = 2    ;
tEnd            = 20    ;
% sIntegration     ='Newmark';
% sIntegration     ='RK4';
sIntegration     ='ode23';
% sIntegration     ='ode45';
% dt              = 0.005 ; % Works with Newmark integration, with or without soil, as long as soil is not in force
dt              = 0.1 ; % Works with Newmark integration, with or without soil, as long as soil is not in force
% dt              = 0.000005         ; % Works with Newmark integration
% dt              = 0.000001       ; % Works for RK4 without CB reduction
% dt              = 0.0005       ; % Works for RK4 with CB reduction

% --- Tube L100-D7-t60   Beam - Structural data
L   = 100      ;
EI0 = 1.654e+12; % [Nm2]
m0  = 1.026e+04; %[kg/m]
A=1;
rho=m0;

% --- Constant for this script
ndofPerNode = 2 ; % Degrees of Freedom per Node: 2 for 2D beams

% vt=linspace(0,20,100);
% Fz=fTowerTopExcitation(vt);
% figure,plot(vt,Fz);
% return
%%

% --------------------------------------------------------------------------------}
%% --- Body 0: Soil element 
% --------------------------------------------------------------------------------{
if bUseSoilElement
    if ~isequal(BC,'free-free')
        error('[FAIL] When using a soil element, the BC should be `free-free`')
    end
    if bCBReduction
        if bSoilAsForce
            error('[FAIL] When using a soil element with CB, the soild cannot be a force (requires nodal displacements and force transfer via transformation matrix')
        end

        if ~isequal(CBLeader,'top-bottom')
            error('[FAIL] When using a soil element with CB, the CB leaders should be `top-bottom`')
        end
    end
    dxFEM = L/(nel);
    fprintf('dxFEM %.3f   dxSoil: %.3f\n',dxFEM, dx_soil);

    L  = L-dx_soil; % NOTE: we make the beam shorter
    nel=nel-1     ; %       and we remove one element

    xe   = [0 dx_soil];
    EIe  = [EI0 EI0];
    me   = [m0 m0];
    nel_1=1;
    gravity=0;Mtop=0;xLumped=[];MLumped=[];bStiffeningMtop=0;bStiffeningSelfWeight=0;
    [Me_soil,Ke_soil]=fBeamMatrices2D_2DOF(xe,EIe,me,nel_1,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight);
else
    Ke_soil=[];
end


% --------------------------------------------------------------------------------}
%% --- Body 1 - FEM
% --------------------------------------------------------------------------------{
fprintf('Building FEM model');
gravity = 9.81;
Mtop    = 0   ;
tic();
[MM,KK,KKg,x]=fBeamMatrices2D_2DOF(L,EI0,m0,nel,gravity,Mtop,[],[],true,true);
toc();
fprintf('Number of DOF of full system:    %d\n',size(KK,1));

if bUseSoilElement
    x=x+dx_soil; % we give offset since the bean was set to a smaller length
end


% --------------------------------------------------------------------------------}
%% --- Boundary conditions - CB Reduction 
% --------------------------------------------------------------------------------{
%% --- Applying BC
if isequal(BC,'clamped-free')
    for i=1:ndofPerNode
        KK(1,:)=[]; KK(:,1) = [];
        MM(1,:)=[]; MM(:,1) = [];
    end
    x_nobc=x(2:end);
elseif isequal(BC,'free-free')
    x_nobc=x;
else
    error('Unsupported boundary condition type')
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

%% --- CB reduction
% NOTE: CB reduction performs a reordering of the DOFs, where the Leader DOFs are placed at the top
if bCBReduction
    nDOF_tot=size(MM,1);
    if isequal(CBLeader, 'top')
        Im = [-ndofPerNode+1:0]+size(MM,1);  % Selecting top node (last node) as master 
        xs = x_nobc(1:end-1);
    elseif isequal(CBLeader, 'top-bottom')
        % Selecting top node as master and bottom node only
        Im = [nDOF_tot-1 nDOF_tot 1 2];
        xs = x_nobc(2:end-1);
    else
        error('Unknown CBLeader option %s',CBLeader);
    end

    [fc,Mr,Kr,~,~,Psi_sq, ~, ~, Psi_sm] = fCraigBamptonReduce(Im,nModesCB,MM,KK);
    Psi_sm(abs(Psi_sm)<1e-8)=0;

    fprintf('fc=%8.3f  (first constrained frequency)\n',fc);
    % [fc,Mr,Kr] = fGuyanReduce(Im,MM,KK);
    fprintf('Number of DOF after reduction:  %d\n',size(Kr,1))
    iTT = 1;  % The tower top DOF dispalcement DOF is at 1;
else
    Mr=MM;
    Kr=KK;
    iTT =size(Mr,1)-1;  % The tower top DOF is at nDOF-1;
end
nDOF_tot=size(Mr,1);


% --------------------------------------------------------------------------------}
%% --- Plotting modes 
% --------------------------------------------------------------------------------{
if bPlotModes
    if CBLeader=='top-bottom'
        warning('Not sure if plotting works for CBLeader = top-bottom')
    end
    % --- Guyan Modes
    U_Guy=Psi_sm(1:ndofPerNode:end,:);
    V_Guy=Psi_sm(2:ndofPerNode:end,:);
    figure
    for i=1:length(Im)
        subplot(1,length(Im),i); hold all; box on
        plot(xs,U_Guy(:,i))
        plot(xs,V_Guy(:,i))
        legend('U','V')
        title(sprintf('Guyan mode %d',i))
    end
    % --- CB modes
    U_CB=Psi_sq(1:ndofPerNode:end,:);
    V_CB=Psi_sq(2:ndofPerNode:end,:);
    figure
    for i=1:floor(nModesCB/2)
        subplot(1,floor(nModesCB/2),i); hold all; box on
        plot(xs,U_CB(:,i))
        plot(xs,V_CB(:,i))
        legend('U','V')
        title(sprintf('CB mode %d',i))
    end
end




% --------------------------------------------------------------------------------}
%% --- Damping 
% --------------------------------------------------------------------------------{
Alpha = 6.0905e-03;
Beta  = 1.4950e-03;
Dr=Alpha*Mr+Beta*Kr;





% --------------------------------------------------------------------------------}
%% --- Inserting soil element 
% --------------------------------------------------------------------------------{
if bUseSoilElement
    if bCBReduction
        if isequal(CBLeader, 'top')
            Isoil=[1 2]; 
        elseif isequal(CBLeader, 'top-bottom')
            Isoil=[3 4]; % nodes 1 2 of original model are nodes 3 4 or reduced model
        end
    else
        Isoil=[1 2]; 
    end
    if ~bSoilAsForce
        % The soil element is inserted in the system stiffness
        Kr(Isoil,Isoil)= Kr(Isoil,Isoil) + Ke_soil(3:4,3:4);
    else
        % The soil element will be included as a force
    end
else
    Isoil=[];
end


% --------------------------------------------------------------------------------}
%% ---  Time Loop
% --------------------------------------------------------------------------------{
vt = 0:dt:tEnd; % NOTE: Purely for storing outputs


%% --- LU decomposition
[M_L,M_U] = lu(Mr);

% ODE
ode_options=struct();
ode_options.OutputFcn=@fodeProgressBar;

nit=0; tic()
if isequal(sIntegration,'Newmark')
    [vt,Y] = fodeNewmark(@(t,x,xp)fTowerTopMDKR(t,x,xp,Mr,Kr,Dr,iTT, @fTowerTopExcitation, Ke_soil,Isoil),vt,zeros(2*nDOF_tot,1));
elseif isequal(sIntegration,'RK4')
    [vt,Y] = fodeRK4(@(t,y)fTowerTopYDot(t,y,M_L,M_U,Kr,Dr,iTT, @fTowerTopExcitation, Ke_soil,Isoil),vt,zeros(2*nDOF_tot,1));
elseif isequal(sIntegration,'ode23')
    [vt,Y] = ode23(@(t,y)fTowerTopYDot(t,y,M_L,M_U,Kr,Dr,iTT, @fTowerTopExcitation, Ke_soil,Isoil),vt,zeros(2*nDOF_tot,1));
elseif isequal(sIntegration,'ode45')
    [vt,Y] = ode45(@(t,y)fTowerTopYDot(t,y,M_L,M_U,Kr,Dr,iTT, @fTowerTopExcitation, Ke_soil,Isoil),vt,zeros(2*nDOF_tot,1));
end
toc();

%% Computing tip-displacement
if bCBReduction
    tip_disp=Y(:,1); % DOF have been rearranged, tip is now at top
else
    tip_disp=Y(:,nDOF_tot-1); % tip is at bottom since order starts from bottom
end

%% Plotting tip-displacement
M_ref = dlmread('RefDisp.csv');
figure, hold all, box on; grid on;
plot(vt   , tip_disp   , '-', 'LineWidth',2)
plot(M_ref(:,1), M_ref(:,2) , 'k.')
xlabel('time [s]'); ylabel('Tower top displacement [m]')
legend('simulation','ref','Location','NorthWest')
title(' Tower top displacement')
xlim([0 max(vt)])
