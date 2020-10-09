function [freq,x,ModesV,ModesK,p] = fUniformBeamTheoryTorsion(Type,G,Kt,Ip,rho,A,L,varargin)
%  returns torsional Modes for a uniform beam 
%
% References:
%  Inman : Engineering variation
% 
%  Author: E. Branlard
%  Date:   06/2017   
%% TEST FUNCTION
if nargin==0
    L   = 100                      ;
    G   = 79.3e9                   ;% Shear modulus. Steel: 79.3  [Pa] [N/m^2]
    D   = 8                        ;
    t   = 0.045                    ;
    A   = pi*( (D/2)^2 - (D/2-t)^2);
    rho = 7850                     ;
    Ip = pi/32*(D^4-(D-2*t)^4); % Polar second moment of area [m^4]
    Kt = pi/64*(D^4-(D-2*t)^4); % Torsion constant
    [freq,x,V,~] = fUniformBeamTheoryTorsion('torsion-unloaded-clamped-free',G,Kt,Ip,rho,A,L,'norm','tip_norm');
    nModesPlot=4;
    close all;
    figure,hold all;box on;
    for i=1:nModesPlot 
        plot(x,V(i,:))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d f=%.8f \n',i,freq(i))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d vmid=%.8f \n',i,V(i,51))
    end

    legds=arrayfun(@(i,f)sprintf('Mode %d f=%4.1f',i,f),1:nModesPlot,freq(1:nModesPlot),'UniformOutput',false);
    legend(legds);
    return
end


%% Optional arguments as Key-values
p = fInputParser();
p.addParameter('x'   ,[],@isnumeric);
p.addParameter('norm','',@ischar   );
p.addParameter('nModes',[],@isnumeric);
p.parse(varargin{:});
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
ModesV = [];
% ModesV = [];
% ModesK = [];


% Dimensionless spanwise position
x0=x/L;

s=strsplit(Type,'-');

if isequal(lower(s{1}),'torsion')
    % --------------------------------------------------------------------------------}
    %% --- PURE TORSION 
    % --------------------------------------------------------------------------------{
    if isequal(lower(s{2}),'unloaded')
        switch(lower(Type))
            case 'torsion-unloaded-clamped-free'
                c=sqrt(G*Kt/(rho*Ip));
                freq   = zeros(1,p.nModes)       ;
                ModesV = NaN(p.nModes,length(x0));
                for j=1:p.nModes
                    omega_j = c/L*(pi/2 + (j-1)*pi);
                    freq(j) = omega_j/(2*pi)       ;
                    ModesV(j,:) = sin(omega_j/c * x);
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


%% Normalization of modes
if isequal(p.norm,'tip_norm')
    for i=1:size(ModesV,1)
        fact=1 / ModesV(i,end);
        ModesV(i,:) = ModesV(i,:)*fact;
    end
end


%% Computation of derivatives if no analytical functions
if exist('fgradient_regular')
    dx=x(2)-x(1);
    ModesK=zeros(size(ModesV));
    for i=1:size(ModesV,1)
        ModesK(i,:)=fgradient_regular(ModesV(i,:),4,dx);
    end
else
    ModesK=[];
end
