function [Mr,Kr,f,x,U,V,K,MM,KK]=fClampedBeamFEM2D(x0,EI0,m0,nel,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight)
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
% OUTPUTS
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

[f,x,U,V,K,MM,KK,Mr,Kr]=fClampedBeamModes2D_2DOF(x0,EI0,m0,nel,gravity,Mtop,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight);
