function [freq,x,ModesU,p] = fUniformBeamTheoryLongi(Type,E,rho,A,L,varargin)
%  returns logitudinals Modes for a uniform beam 
%
% References:
% 
%  Author: E. Branlard
%  Date:   06/2017   
%% TEST FUNCTION
if nargin==0
    L   = 100                      ;
    D   = 8                        ;
    t   = 0.045                    ;
    A   = pi*( (D/2)^2 - (D/2-t)^2);
    E   = 210e9                    ;% Young modulus [Pa] [N/m^2]
    rho = 7850                     ;
    [freq,x,U,p] = fUniformBeamTheoryLongi('longi-unloaded-clamped-free',E,rho,A,L,'norm','tip_norm');
    nModesPlot=4;
    close all;
    figure,hold all;box on;
    for i=1:nModesPlot 
        plot(x,U(i,:))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d f=%.8f \n',i,freq(i))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d umid=%.8f \n',i,U(i,51))
    end


    legds=arrayfun(@(i,f)sprintf('Mode %d f=%4.1f',i,f),1:nModesPlot,freq(1:nModesPlot),'UniformOutput',false);
    legend(legds);
    return
end


%% Optional arguments as Key-values
p = inputParser;
% --- Key, value parameters
addParameter(p,'x'   ,[],@isnumeric)
addParameter(p,'norm',[],@ischar   )
addParameter(p,'nModes',[],@isnumeric)
parse(p,varargin{:});
p=p.Results;
%
if isempty(p.x)
    x= linspace(0,L,101);
else
    x=p.x;
    if max(x)~=L; error('Max of x should be equal to L'); end;
end;
if isempty(p.nModes)
    p.nModes=20;
end

%
freq   = [];
ModesU = [];
% ModesV = [];
% ModesK = [];


% Dimensionless spanwise position
x0=x/L;

s=strsplit(Type,'-');

if isequal(lower(s{1}),'longi')
    % --------------------------------------------------------------------------------}
    %% --- PURE TORSION 
    % --------------------------------------------------------------------------------{
    if isequal(lower(s{2}),'unloaded')
        switch(lower(Type))
            case 'longi-unloaded-clamped-free'
                m=rho*A; % [kg/m]
                c=sqrt(E*A/m);
                freq   = zeros(1,p.nModes)       ;
                ModesU = NaN(p.nModes,length(x0));
                for j=1:p.nModes
                    omega_j = c/L*(pi/2 + (j-1)*pi);
                    freq(j) = omega_j/(2*pi)       ;
                    ModesU(j,:) = sin(omega_j/c * x);
                end

            otherwise
                error('unknown type %s',Type)
        end % switch
    else
        error('Unknown %s',Type);
    end % transverse loaded
else
    error('Unknown %s',Type);
end

%% Computation of derivatives if no analytical functions
% V=fgradient_regular(U(i,:),4,dx);
% K=fgradient_regular(V(i,:),4,dx);

%% Going back to physical dimension
% x=x0*L;
% ModesV = ModesV/L;
% ModesK = ModesK/L^2;

%% Normalization of modes
if isequal(p.norm,'tip_norm')
    for i=1:size(ModesU,1)
        fact=1 / ModesU(i,end);
        ModesU(i,:) = ModesU(i,:)*fact;
    end
end


