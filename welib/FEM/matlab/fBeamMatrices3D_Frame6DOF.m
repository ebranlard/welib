function [MM,KK,x]=fBeamMatrices3D_Frame6DOF(x0,EIx0,EIy0,EIz0,EA0,m0,A0,E,G,Kt,nel)

% Returns the mass and stiffness FEM matrix of a beam represented with nel Frame elements 
%
% NOTE: input values x0, EI0 m0 can be vectors or scalars.
%       If they are scalars, then a beam with constant properties and of length L=x0 is used;
%       If they are vectors, then linear interpolation is used. The dimension of the inputs does not need to match nel
% 
% INPUTS
%   x0   : (1xn) Span vector of the beam [m]
%   EI0  : (1xn) Elastic Modulus times Second Moment of Area of cross section [Nm2]
%   m0   : (1xn) Mass per length along the beam [kg/m]
%   nel  : Number of elements
%
% OUTPUTS
%   MM: (nDOF x nDOF)  Mass matrix
%   KK: (nDOF x nDOF)  Stiffness matrix
%   x : (1 x nel)   Span vector
%
% AUTHOR: E. Branlard 
%

% --------------------------------------------------------------------------------}
%% --- Optional Input  
% --------------------------------------------------------------------------------{
% TODO:
% if ~exist('bStiffeningMtop'    ,'var');    bStiffeningMtop=true; end
% if ~exist('bStiffeningSelfWeight','var');  bStiffeningSelfWeight=true; end;
% if ~exist('gravity','var');  gravity=9.81; end;
% if ~exist('Mtop','var');     Mtop=0; end;
% if ~exist('MLumped','var');  MLumped=[]; end;
% if ~exist('xLumped','var');  xLumped=[]; end;


if length(x0)==1
    % Constant beam properties
    x0   = [0 x0]    ;
    EIx0 = [1 1]*EIx0;
    EIy0 = [1 1]*EIy0;
    EIz0 = [1 1]*EIz0;
    EA0  = [1 1]*EA0 ;
    A0   = [1 1]*A0  ;
    m0   = [1 1]*m0  ;
end

%% Parameters
% Beam Element:
ndof=6; % Degrees of Freedom per Node
nnel=2;          % number of nodes per element
%
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

% --------------------------------------------------------------------------------}
%% --- Mesh and Geometry 
% --------------------------------------------------------------------------------{
P=zeros(3,nel+1);
P(1,:)=linspace(min(x0),max(x0),nel+1);
x=P(1,:);

% --------------------------------------------------------------------------------}
%% --- Computing element properties
% --------------------------------------------------------------------------------{
% Calculating where X is the position of element from Fixed End
Mass = zeros(1,nel);
Ae   = zeros(1,nel);
me   = zeros(1,nel);
Le   = zeros(1,nel);
EAe  = zeros(1,nel);
EIxe = zeros(1,nel);
EIye = zeros(1,nel);
EIze = zeros(1,nel);
Kte  = zeros(1,nel);

for i=1:nel
    % Element Points
    P1   = x(i);
    P2   = x(i+1);
    Pmid = (P2+P1)/2;
    % Length of Element
    Le(i)=norm(P2-P1); 
    % Structural data of element interpolated from inputs
    EIxe(i)  = interp1(x0,EIx0,Pmid(1)) ;
    EIye (i) = interp1(x0,EIy0 ,Pmid(1));
    EIze (i) = interp1(x0,EIz0 ,Pmid(1));
    EAe  (i) = interp1(x0,EA0  ,Pmid(1));
    me   (i) = interp1(x0,m0   ,Pmid(1));
    Ae   (i) = interp1(x0,A0   ,Pmid(1));
    Mass(i)  = me(i)*Le(i)              ; % TODO add concentrated mass
    Kte(i)   = Kt                       ;
end
% --------------------------------------------------------------------------------}
%% --- Building 
% --------------------------------------------------------------------------------{
MM =  zeros(sdof,sdof);
KK =  zeros(sdof,sdof);
% Loop on elements
for i=1:nel
    DOFindex=fElementDOFIndex(i,nnel,ndof); % 1 x ndof*nnel
    % --- Element matrix
    [ke,me]=fElementMatricesFrame3D_6DOF(E,G,Kte(i),EAe(i),EIxe(i),EIye(i),EIze(i),Le(i),Ae(i),Mass(i));
    % TODO: geometric stiffness
    % --- Build global matrices
    MM =fBuildGlobalMatrix(MM ,me  ,DOFindex);
    KK =fBuildGlobalMatrix(KK ,ke  ,DOFindex);
end

