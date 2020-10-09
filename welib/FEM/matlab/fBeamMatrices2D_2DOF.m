function [MM,KK,KKg,x]=fBeamMatrices2D_2DOF(x0,EI0,m0,nel,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight);
% Returns the mass and stiffness FEM matrix of a beam represented with nel elements 
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
%   KKg:(nDOF x nDOF)  Geometric stiffness matrix
%   x : (1 x nel)   Span vector
%
% AUTHOR: E. Branlard
%
%% Test Function
if nargin==0
    nel=100; % Number of Elements
    % --- Tube Beam - Structural data
    L  = 100                  ;
    EI = 1.868211939147334e+12;
    m  = 8.828201296825122e+03;
    [MM,KK,KKg,x]=fBeamMatrices2D_2DOF(L,EI,m,nel);
    return
end

% --------------------------------------------------------------------------------}
%% --- Optional Input  
% --------------------------------------------------------------------------------{
if ~exist('bStiffeningMtop'    ,'var');    bStiffeningMtop=true; end
if ~exist('bStiffeningSelfWeight','var');  bStiffeningSelfWeight=true; end;
if ~exist('gravity','var');  gravity=9.81; end;
if ~exist('Mtop','var');     Mtop=0; end;
if ~exist('MLumped','var');  MLumped=[]; end;
if ~exist('xLumped','var');  xLumped=[]; end;

if length(xLumped)>0
    error('So far this script is only for Mtop not for Mlump')
end

if length(x0)==1
    % Constant beam properties
    x0  = [0 x0];
    EI0 = [1 1]*EI0;
    m0  = [1 1]*m0;
end

% --------------------------------------------------------------------------------}
%% --- Mesh and Geometry 
% --------------------------------------------------------------------------------{

if isempty(nel)
    nel=length(x0)-1;
    P=zeros(3,nel+1);
    P(1,:)=x0;
    x=P(1,:);
else
    % Using a linear mesh
    P=zeros(3,nel+1);
    P(1,:)=linspace(min(x0),max(x0),nel+1);
    x=P(1,:);
end


%% --- Parameters 
% Beam Element:
ndof=2; % Degrees of Freedom per Node
nnel=2; % number of nodes per element
%
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  


% --------------------------------------------------------------------------------}
%% --- Computing element properties, and axial force 
% --------------------------------------------------------------------------------{
Mass = zeros(1,nel);
EI   = zeros(1,nel);
m    = zeros(1,nel);
Le   = zeros(1,nel);
for i=1:nel
    % Element Points
    P1   = x(i);
    P2   = x(i+1);
    Pmid = (P2+P1)/2;
    % Length of Element
    Le(i)=norm(P2-P1); 
    % Structural data of element interpolated from inputs
    EI(i)  = interp1(x0,EI0 ,Pmid(1));
    m (i)  = interp1(x0,m0,Pmid(1));
    Mass(i)=m(i)*Le(i); % TODO add concentrated mass
end
% Self Weigth contribution to Pacc
if ~exist('flip')
    Pacc_SW = - gravity*(fliplr(cumsum(Mass))); % NOTE: negative for compression
else
    Pacc_SW = - gravity*(flip(cumsum(Mass))); % NOTE: negative for compression
end
Pacc_MT = - gravity*Mtop*ones(1,nel)    ; % NOTE: negative for compression
Pacc    = zeros(1,nel)                  ;

% TopMass contribution to Pacc
if bStiffeningMtop
    Pacc=Pacc+Pacc_MT;
end
if bStiffeningSelfWeight
    Pacc=Pacc+Pacc_SW;
end



% --------------------------------------------------------------------------------}
%% --- Building 
% --------------------------------------------------------------------------------{
MM =  zeros(sdof,sdof);
KK =  zeros(sdof,sdof);
KKg=  zeros(sdof,sdof);
% Loop on elements
for iel=1:nel
    DOFindex=fElementDOFIndex(iel,nnel,ndof); % 1 x ndof*nnel
    M_e = m(iel)*Le(iel); % Element mass in kg:  M_e=rho A L = m L
    % --- Element matrix
    [ke1,me]  = fElementMatricesBeam2D_2DOF(EI(iel),Le(iel),M_e,1)  ;
    [kg]      = fElementKNormalForce2D_2DOF(Le(iel),Pacc(iel))          ;
    % --- Build global matrices
    MM =fBuildGlobalMatrix(MM ,me  ,DOFindex);
    KK =fBuildGlobalMatrix(KK ,ke1 ,DOFindex);
    KKg=fBuildGlobalMatrix(KKg,kg  ,DOFindex);
end
%% Adding concentrated masses
MM(end-1,end-1)=MM(end-1,end-1)+Mtop;
