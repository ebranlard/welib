function [f,x,U,V,K]=fClampedBeamRotatingModes(x0,EI0,m0,Omega,nel);
% Use FEM to get the modes of a rotating or non-rotating cantilever beam
%
% NOTE: input values can be vectors or scalars.
%       If they are scalars, then a beam with constant properties and of length L=x0 is used;
%       If they are vectors, then linear interpolation is used. The dimension of the inputs does not need to match nel
% 
% INPUTS
%   x0   : (1xn) Span vector of the beam [m]
%   EI0  : (1xn) Elastic Modulus times Second Moment of Area of cross section [Nm2]
%   m0   : (1xn) Mass per length along the beam [kg/m]
%   Omega: Rotational speed of the beam
%   nel  : Number of elements
%
% OUTPUS
%   f : (1 x nMl)   Modes Frequencies
%   x : (1 x nel)   Span vector
%   U : (nM x nel)  Modes Shapes Displacements
%   V : (nM x nel)  Modes Shapes Slopes
%   K : (nM x nel)  Modes Shapes Curvature
%
%


%% Test Function
if nargin==0
    nel=100; % Number of Elements
    Omega =6; % [rad/s]
    % --- Kane Beam - Structural data
    x0   = 10;
    EI0  = 14000;
    m0   = 1.2;
    [f,x,U,V,K]=fClampedBeamRotatingModes(x0,EI0,m0,Omega,nel);
    x_th=x;
    [f_th,~,U_th,V_th,K_th] = fUniformBeamTheory('transverse-unloaded-clamped-free',EI0(1),m0(1),1,max(x0),'x',x_th,'norm','tip_norm');
    %% Comparison of Frequencies
    for i=1:min(8,nel)
         fprintf('f%d=%8.3f  -   f=%8.3f   -   df=%9.5f\n',i,f(i),f_th(i),f_th(i)-f(i))
    end
    %% Plotting Mode Shapes
    vColrs=lines;
    figure(1),clf,hold all,box on;
    for i=1:min(4,nel)
        plot(x_th,U_th(i,:),'-'  ,'Color',vColrs(i,:));
        plot(x   ,U   (i,:),'.--','Color',vColrs(i,:));
    end
    legend('Non rotating','Rotating'); title('Mode shapes - U')
    %% Plotting Mode Shapes Slopes
    figure(2),clf,hold all,box on;
    for i=1:min(4,nel)
        plot(x_th,V_th(i,:),'-'  ,'Color',vColrs(i,:));
        plot(x   ,V(i,:)   ,'.--','Color',vColrs(i,:));
    end
    legend('Non rotating','Rotating'); title('Mode shapes slopes - V')
    %% Plotting Mode Shapes Curvatures
    figure(3),clf,hold all,box on;
    for i=1:min(4,nel)
        plot(x_th,K_th(i,:),'-' ,'Color',vColrs(i,:));
        plot(x   ,K(i,:)  ,'.--','Color',vColrs(i,:));
    end
    legend('Non rotating','Rotating'); title('Mode shapes curvatures - K')
end

% --------------------------------------------------------------------------------}
%% --- Input  
% --------------------------------------------------------------------------------{
if length(x0)==1
    % Constant beam properties
    x0  = [0 x0];
    EI0 = [1 1]*EI0;
    m0  = [1 1]*m0;
end





%% --- Parameters 
% Beam Element:
ndof=2; % Degrees of Freedom per Node
nnel=2; % number of nodes per element
%
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

% --------------------------------------------------------------------------------}
%% --- Mesh and Geometry 
% --------------------------------------------------------------------------------{
P=zeros(3,nel+1);
P(1,:)=linspace(min(x0),max(x0),nel+1);
x=P(1,:);
% Calculating A value for each element
AA=zeros(1,nel);
for i=1:nel
    for j=i:nel
        AA(i)= AA(i)+((x(j+1))^2)-((x(j))^2);
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
    node1=iel;      % starting node number for element 'iel'
    node2=iel+1;    % ending node number for element 'iel'
    % Element Points
    P1   = P(:,node1) ;
    P2   = P(:,node2) ;
    Pmid = (P2+P1)/2;
    % Length of Element
    Le=norm(P2-P1); 
    % Structural data of element interpolated from inputs
    EI  = interp1(x0,EI0,Pmid(1));
    m   = interp1(x0,m0 ,Pmid(1));
    M_e = m*Le; % Element mass in kg:  M_e=rho A L = m L
    % --- Element matrix
    [ke1,me]      = fElementMatricesBeam2DOF(EI,Le,M_e,1)  ;
    [ke2,ke3,ke4] = fElementKeRotatingBeam2DOF(Le);
    ke=ke1+ (m*(Omega^2)*AA(iel)/2)*ke2 -(m*Omega^2*x(iel)/2)*ke3 -(m*Omega^2/2)*ke4;
    % --- Build global matrices
    MM =fBuildGlobalMatrix(MM ,me  ,DOFindex);
    KK =fBuildGlobalMatrix(KK ,ke  ,DOFindex);
end

%% --- Apply BC
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
f    = sqrt(Omega2)/(2*pi);
% Freq_NonDim_ = sqrt(Omega2)/sqrt((E*I)/(m * (L^4)));
% Omega = lambda * sqrt((E*I)/(m * (L^4)));% Rotating Speed in radian/sec
% lambda =Omega / sqrt((E*I)/(m * (L^4)));% Rotating Speed in radian/sec

%% Mode shapes - Normalized to unit tip
U=Q(1:ndof:end,:)';
V=Q(2:ndof:end,:)';
% Adding 0s that were removed due to BC
U=[zeros(size(U,1),1) U];
V=[zeros(size(V,1),1) V];
K = nan(size(U));

dx=x(2)-x(1);
for i=1:size(U,1)
    fact=1 / U(i,end);
    U(i,:) = U(i,:)*fact;
    V(i,:) = V(i,:)*fact;
    % Computing K
    K(i,:) = fgradient_regular(V(i,:),4,dx);
end
