function [f,x,U,V,K,MM,KK,Mr,Kr]=fClampedBeamModes2D_2DOF(x0,EI0,m0,nel,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight)
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
    % --- Tube Beam - Structural data
    L  = 100                  ;
    EI = 1.868211939147334e+12;
    m  = 8.828201296825122e+03;
    Mtop=1*m*L;
    [f,x,U,V,K]=fClampedBeamModes2D_2DOF(L,EI,m,nel,9.81,Mtop,[],[],1,1);
    x_th=x;
    addpath('../BeamTheory/');
    [f_th,~,~,~,~] = fUniformBeamTheory('transverse-unloaded-topmass-clamped-free',EI,m,1,L,'x',x_th,'norm','tip_norm','Mtop',Mtop);
    %% Comparison of Frequencies
    for i=1:min(8,nel)
         fprintf('f%d=%8.3f  -   f=%8.3f   -   df=%9.5f\n',i,f(i),f_th(i),f_th(i)-f(i))
    end
    return
end

% --------------------------------------------------------------------------------}
%% --- Input  
% --------------------------------------------------------------------------------{
if ~exist('bStiffeningMtop'    ,'var');    bStiffeningMtop=true; end
if ~exist('bStiffeningSelfWeight','var');  bStiffeningSelfWeight=true; end;
if ~exist('gravity','var');  gravity=9.81; end;
if ~exist('Mtop','var');     Mtop=0; end;
if ~exist('MLumped','var');  MLumped=[]; end;
if ~exist('xLumped','var');  xLumped=[]; end;

if ~isempty(xLumped)
    error('So far this script is only for Mtop, Mlumped to be implemented')
end

if length(x0)==1
    % Constant beam properties
    x0  = [0 x0];
    EI0 = [1 1]*EI0;
    m0  = [1 1]*m0;
end





%% --- Parameters 
% Beam Element:
ndof=2; % Degrees of Freedom per Node
%
%% --- Building mass and stiffness matrix
[MM,KK,KKg,x]=fBeamMatrices2D_2DOF(x0,EI0,m0,nel,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight);
KK=KK+KKg;

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
[Omega2,Isort]=sort(Omega2);
Q=Q(:,Isort);
f= sqrt(Omega2)/(2*pi);

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
nSpan=size(U,2);

dx=x(2)-x(1);
for i=1:size(U,1)
    fact=1 / U(i,end);
    U(i,:) = U(i,:)*fact;
    V(i,:) = V(i,:)*fact;
    % Computing K
    if nSpan>=5
        K(i,:) = fgradient_regular(V(i,:),4,dx);
    elseif nSpan>=3
        K(i,:) = fgradient_regular(V(i,:),2,dx);
    else
        K(i,:) = 0*V(i,:);% TODO
    end
end
