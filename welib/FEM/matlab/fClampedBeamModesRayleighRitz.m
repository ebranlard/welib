function [f,x,U,V,K,MM,KK]=fClampedBeamModesRayleighRitz(x0,EI0,m0,nel,nModes,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight);
% Use RayleighRitz to get the modes of a rotating or non-rotating cantilever beam
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
    nModes=10; % Number of Modes
    % --- Tube Beam - Structural data
    L  = 100                  ;
    EI = 1.868211939147334e+12;
    m  = 8.828201296825122e+03;
    Mtop=0*m*L;
    [f,x,U,V,K]=fClampedBeamModesRayleighRitz(L,EI,m,nel,nModes,9.81,Mtop,[],[],0,0);
    x_th=x;
    [f_th,~,U_th,V_th,K_th] = fUniformBeamTheory('transverse-unloaded-topmass-clamped-free',EI,m,1,L,'x',x_th,'norm','tip_norm','Mtop',Mtop);

    %% Comparison of Frequencies
    for i=1:min(4,nModes)
         fprintf('f%d=%8.4f  -   f=%8.4f   -   df=%9.5f\n',i,f(i),f_th(i),f_th(i)-f(i))
    end
    figure,hold all
    for i=1:min(4,nModes)
        plot(x_th,U_th(i,:),'k.')
        plot(x,U(i,:))
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

if ~isempty(xLumped); error('So far this script is only for Mtop, lumped masses not implemented'); end
if ~isempty(MLumped); error('So far this script is only for Mtop, lumped masses not implemented'); end

if length(x0)~=1
    error('Rayleigh Ritz method for now works with constant properties. (Otherwise, use FEM to get mode shapes..)');
end
% Constant beam properties
x0  = [0 x0];
EI0 = [1 1]*EI0;
m0  = [1 1]*m0;




% --------------------------------------------------------------------------------}
%% --- Mesh and Geometry 
% --------------------------------------------------------------------------------{
P=zeros(3,nel+1);
P(1,:)=linspace(min(x0),max(x0),nel+1);
x=P(1,:);
m  = interp1(x0,m0,x);
EI = interp1(x0,EI0,x);

% --------------------------------------------------------------------------------}
%% --- Axial force 
% --------------------------------------------------------------------------------{
Pacc_SW = fcumtrapzlr(x, -m * gravity) ;
Pacc_MT = -Mtop * gravity*ones(size(x));
Pacc    = zeros(1,nel+1)               ;

% TopMass contribution to Pacc
if bStiffeningMtop
    Pacc=Pacc+Pacc_MT;
end
if bStiffeningSelfWeight
    Pacc=Pacc+Pacc_SW;
end

% --------------------------------------------------------------------------------}
%% --- Mode definition 
% --------------------------------------------------------------------------------{
x=linspace(0,max(x0),nel+1);
[f_th,x_th,U_th,V_th,K_th] = fUniformBeamTheory('transverse-unloaded-topmass-clamped-free' ,EI0(1),m0(1),1,max(x0),'x',x,'norm','tip_norm','Mtop',Mtop);
%
U0=U_th(1:nModes,:);
V0=V_th(1:nModes,:);
K0=K_th(1:nModes,:);


% --------------------------------------------------------------------------------}
%% --- Building 
% --------------------------------------------------------------------------------{
MM0=  zeros(nModes,nModes);
KK0=  zeros(nModes,nModes);
KKg=  zeros(nModes,nModes);
%% Building matrices
for i=1:nModes
    MM0(i,i)  = trapz(x, m   .* U0(i,:).*U0(i,:));
    KK0(i,i)  = trapz(x, EI  .* K0(i,:).*K0(i,:));
    for j=1:nModes
        KKg(i,j) = trapz(x, Pacc .* V0(i,:).*V0(j,:));
    end
end
KK=KK0+KKg;
%% Adding concentrated masses
MM=MM0+eye(nModes)*Mtop;

%% --- Solve EVA
[Q,Lambda]=eig(KK,MM);
Omega2=diag(Lambda);
[Omega2,Isort]=sort(Omega2);
Q=Q(:,Isort);
f= sqrt(Omega2)/(2*pi);

%% Mode shapes - Normalized to unit tip
U=zeros(nModes,nel+1);
V=zeros(nModes,nel+1);
K=zeros(nModes,nel+1);
for i=1:nModes
    for j=1:nModes
        U(i,:)=U(i,:)+ Q(j,i)*U0(j,:);
        V(i,:)=V(i,:)+ Q(j,i)*V0(j,:);
        K(i,:)=K(i,:)+ Q(j,i)*K0(j,:);
    end
    fact=1/U(i,end);
    U(i,:)=U(i,:)*fact;
    V(i,:)=V(i,:)*fact;
    K(i,:)=K(i,:)*fact;
end



