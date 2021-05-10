%% Documentation
%FEA to calculate Rotating and Non-Rotating Natural Frequencies of
% cantilevered uniform beam
%% Init
clear all; clc; close all;
addpath('../BeamTheory')

%% Parametesr
nel=10; % Number of Elements
% 
Omega =0; % [rad/s] <<<<<<<<<<<<<<<<<<<<<<<<<<<

% --- Kane Beam - Structural data
% L = 10  ; % Radius of Blade in m
% m = 1.2 ; % Mass per Unit Length of Blade in kg/m
% E = 7e10;
% I = 2e-7;
% A=1;
% rho=m/A;

% --- Ganguli Beam - Structual Data
% L = 10  ; % Radius of Blade in m
% m = 5 ; % Mass per Unit Length of Blade in kg/m
% E = 100000;
% I = 1;
% A=1;
% rho=m/A;

% --- Tube
L   = 100                      ;
E   = 210e9                    ;
D   = 8                        ;
t   = 0.045                    ;
A   = pi*( (D/2)^2 - (D/2-t)^2);
I   = pi/64*(D^4-(D-2*t)^4)  ; % (Nm^2) I_annulus = pi/4 (r2^4 - r1^4) = pi/64 (D2^4-D1^4)
rho = 7850                     ;
m   = rho*A                    ;


% Structural data
x0     =[0,L]; % Span wise position of beam;
m0     = m  * ones(size(x0)); % [kg/m]
EI0    = E*I* ones(size(x0)); % [Nm2]:  Elastic Modulus [N/m2] Times moment of inertia of cross-section [m4]


%% Non rotating theory
x_th=linspace(0,L,nel+1);
[f_th,x_th,U_th,V_th,K_th] = fUniformBeamTheory('transverse-unloaded-clamped-free',E*I,rho,A,L,'x',x_th,'norm','tip_norm');

%% Parameters
% --- Derived parameters
% Omega = lambda * sqrt((E*I)/(m * (L^4)));% Rotating Speed in radian/sec
% lambda =Omega / sqrt((E*I)/(m * (L^4)));% Rotating Speed in radian/sec


% Beam Element:
ndof=2; % Degrees of Freedom per Node
nnel=2;          % number of nodes per element
%
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  



% --------------------------------------------------------------------------------}
%% --- Boundary Conditions (BC)
% --------------------------------------------------------------------------------{
% Clamped-Free BC
% bcdof(1)=1; bcval(1)=0;  % first element deflection u=0
% bcdof(2)=2; bcval(2)=0;  % first element slope v=0

% --------------------------------------------------------------------------------}
%% --- Mesh and Geometry 
% --------------------------------------------------------------------------------{
l=(L/nel); % Length of Each Element in m
% Calculating where X is the position of element from Fixed End
X=0:l:nel*l;
% Calculating A value for each element
AA=zeros(1,nel);
for i=1:nel
    for j=i:nel
        AA(i)= AA(i)+((X(j+1))^2)-((X(j))^2);
    end
end

% --------------------------------------------------------------------------------}
%% --- Building 
% --------------------------------------------------------------------------------{
MM =  zeros(sdof,sdof);
KK =  zeros(sdof,sdof);
% Loop on elements
for iel=1:nel
    DOFindex=fElementDOFIndex(iel,nnel,ndof); % 1 x ndof*nnel

    % --- Element matrix
    [ke1,me]      = fElementMatricesBeam2DOF(E*I,l,l*A*rho,1)  ;
    [ke2,ke3,ke4] = fElementKeRotatingBeam2DOF(l);
    ke=ke1+ (m*(Omega^2)*AA(iel)/2)*ke2 -(m*Omega^2*X(iel)/2)*ke3 -(m*Omega^2/2)*ke4;
    % --- Build global matrices
    MM =fBuildGlobalMatrix(MM ,me  ,DOFindex);
    KK =fBuildGlobalMatrix(KK ,ke  ,DOFindex);
end

% --- Apply BC
% Clamped-Free:  - Removing two first column and row of M and K
Kr=KK;
Mr=MM;
Kr(1,:)=[]; Kr(:,1) = [];
Kr(1,:)=[]; Kr(:,1) = [];
Mr(1,:)=[]; Mr(:,1) = [];
Mr(1,:)=[]; Mr(:,1) = [];

%% --- Solve EVA
[Q,Lambda]=eig(Kr,Mr);
Omega2=diag(Lambda);
bPos=(Omega2>0);
Omega2 = Omega2(bPos);
Q      = Q(:,bPos)   ;
Freq   = sqrt(Omega2)/(2*pi);
fprintf('%d non-positive frequencies discarded \n',sum(~bPos));




%% Comparison of Frequencies
for i=1:min(8,length(Freq));
    fprintf('f%d=%8.3f  -   f=%8.3f   -   df=%9.5f\n',i,Freq(i),f_th(i),f_th(i)-Freq(i))
end
% Natural Frequencies of rotating blade in non-dimensional form
Non_dim_nat_freq = sqrt(Omega2)/sqrt((E*I)/(m * (L^4)));


%% Mode shapes - Normalized to unit tip
U=Q(1:ndof:end,:);
V=Q(2:ndof:end,:);
for i=1:size(U,1)
    fact=1 / U(end,i);
    U(:,i) = U(:,i)*fact;
    V(:,i) = V(:,i)*fact;
end

%th=x_th*L;
%% Plotting Mode Shapes
vColrs=lines;
figure(1),clf,hold all,box on;
for i=1:4
    plot(x_th      ,U_th(i,:),'-','Color',vColrs(i,:));
    plot(X,[0; U(:,i)],'.--','Color',vColrs(i,:));
end
legend('Non rotating','Rotating')


%% Plotting Mode Shapes Slopes
vColrs=lines;
figure(2),clf,hold all,box on;
for i=1:4
    plot(x_th      ,V_th(i,:),'-','Color',vColrs(i,:));
    plot(X,[0; V(:,i)],'.--','Color',vColrs(i,:));
end
legend('Non rotating','Rotating')
