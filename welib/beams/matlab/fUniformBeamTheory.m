function [freq,x,ModesU,ModesV,ModesK,p] = fUniformBeamTheory(Type,EI,rho,A,L,varargin)
%  returns Modes for a uniform beam 
%
% References:
%  Inman : Engineering variation
% 
%  Author: E. Branlard
%  Date:   06/2017   
%% TEST FUNCTION
if nargin==0
    L  = 100                  ;
    EI = 1.868211939147334e+12;
    m  = 8.828201296825122e+03;
    [freq,x,U,V,K,~,~] = fUniformBeamTheory('transverse-unloaded-clamped-free',EI,m,1,L,'norm','tip_norm','x',x);
    nModesPlot=6;
    close all;
    figure,hold all;box on;
    for i=1:nModesPlot 
        fprintf('Mode %d f=%.8f \n',i,freq(i))
        plot(x,U(i,:))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d umid=%.8f \n',i,U(i,51))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d vmid=%.8f \n',i,V(i,51))
    end
    for i=1:nModesPlot 
        fprintf('Mode %d kmid=%.8f \n',i,K(i,51))
    end
    legds=arrayfun(@(i,f)sprintf('Mode %d f=%4.1f',i,f),1:nModesPlot,freq(1:nModesPlot),'UniformOutput',false);
    legend(legds);
    return
end


%% Optional arguments as Key-values
p =fInputParser();
% --- Key, value parameters
addParameter(p,'w'   ,[],@isnumeric)
addParameter(p,'x'   ,[],@isnumeric)
addParameter(p,'Mtop',[],@isnumeric)
addParameter(p,'norm',[],@ischar   )
parse(p,varargin{:});
p=p.Results;
%
if isempty(p.x)
    x= linspace(0,L,101);
else
    x=p.x;
    if max(x)~=L; error('Max of x should be equal to L'); end;
end;

%
freq   = [];
ModesU = [];
ModesV = [];
ModesK = [];


% Dimensionless spanwise position
x0=x/L;

s=strsplit(Type,'-');

if isequal(lower(s{1}),'transverse')
    if isequal(lower(s{2}),'unloaded')
        % --- "Theory" (clamped-free, vertical, no weight)
        % See Inman, p.335 or Nielsen1 p129
        switch(lower(Type))
            case 'transverse-unloaded-clamped-free'
                % NOTE: cosh(beta_n)cos(beta_n) =-1
                %    sigma_n = [ sinh(beta_n) - sin(beta_n) ]/[cosh(beta_n) + cos(beta_n)]
                %    for j>5, a good approx is B(j) = (2*j-1)pi/2  and S(j)=1;
                %B  = [1.87510407, 4.69409113, 7.85475744,10.99554073,14.13716839, (2*6-1)*pi/2];
                %S  = [0.734095514 1.018467319 0.999224497 1.000033553 0.999998550 1];
                B=zeros(1,12); % Beta
                for i=1:length(B)
                    B(i)=fzero(@(x) 1+cosh(x)*cos(x),(2*i-1)*pi/2, optimset('TolFun',1e-21));
                end
            case 'transverse-unloaded-topmass-clamped-free'
                % The geometrical stiffning is not accounted for here
                if isempty(p.Mtop); error('Please specify value for Mtop for %s',Type); end;
                Mtop = p.Mtop    ;
                M    = rho *A * L;
                B=zeros(1,12); % Beta
                for i=1:length(B)
                    B(i)=fzero(@(x) 1+cosh(x)*cos(x) - x*Mtop/M *(sin(x)*cosh(x)-cos(x)*sinh(x)),(2*i-1)*pi/2, optimset('TolFun',1e-21));
                end
            otherwise
                error('unknown type %s',Type)
        end % swtich
        %S  = ( sinh(B)-sin(B) ) ./ ( cosh(B) + cos(B));  % Sigma
        %C  = ( cosh(B)+cos(B) ) ./ ( sinh(B) + sin(B));  % Sigma
        SS  = sinh(B)+sin(B);
        CC  = cosh(B)+cos(B); 

        % Frequency
        freq = (B/L).^2/(2*pi)*sqrt(EI /(rho*A)); 

        % --- Mode shapes
        ModesU=zeros(length(B),length(x0));
        for i=1:length(B)
            ModesU(i,:) =         SS(i)*(cosh(B(i)*x0)-cos(B(i)*x0)) - CC(i)*(sinh(B(i)*x0)-sin(B(i)*x0)) ;
            ModesV(i,:) = B(i)  *(SS(i)*(sinh(B(i)*x0)+sin(B(i)*x0)) - CC(i)*(cosh(B(i)*x0)-cos(B(i)*x0)));
            ModesK(i,:) = B(i)^2*(SS(i)*(cosh(B(i)*x0)+cos(B(i)*x0)) - CC(i)*(sinh(B(i)*x0)+sin(B(i)*x0)));
%             ModesU(i,:)  =        cosh(B(i)*x0)-cos(B(i)*x0) - S(i)*(sinh(B(i)*x0)-sin(B(i)*x0)) ;
%             ModesV(i,:) = B(i)  *(sinh(B(i)*x0)+sin(B(i)*x0) - S(i)*(cosh(B(i)*x0)-cos(B(i)*x0)));
%             ModesK(i,:) = B(i)^2*(cosh(B(i)*x0)+cos(B(i)*x0) - S(i)*(sinh(B(i)*x0)+sin(B(i)*x0)));
%             ModesU(i,:)  =        cosh(B(i)*x0)-cos(B(i)*x0) - C(i)*(sinh(B(i)*x0)-sin(B(i)*x0)) ;
%             ModesV(i,:) = B(i)  *(sinh(B(i)*x0)+sin(B(i)*x0) - C(i)*(cosh(B(i)*x0)-cos(B(i)*x0)));
%             ModesK(i,:) = B(i)^2*(cosh(B(i)*x0)+cos(B(i)*x0) - C(i)*(sinh(B(i)*x0)+sin(B(i)*x0)));
        end
    elseif isequal(lower(s{2}),'loaded')
        switch(lower(Type))
            case 'transverse-loaded-clamped-free'
                if isempty(p.w); p.w= A*rho; end;
                if isempty(p.L); error('Please specify value for L for %s',Type); end;
                B =[1.875,4.694];
                freq = (B/L).^2/(2*pi)*sqrt(EI /w);

            otherwise
                error('unknown type %s',Type)
        end % switch
    else
        error('Unknown %s',Type);
    end % transverse loaded
elseif isequal(lower(s{1}),'torsion')
    % --------------------------------------------------------------------------------}
    %% --- PURE TORSION 
    % --------------------------------------------------------------------------------{
    if isequal(lower(s{2}),'unloaded')
        switch(lower(Type))
            case 'torsion-unloaded-clamped-free'

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
x=x0*L;
ModesV = ModesV/L;
ModesK = ModesK/L^2;

%% Normalization of modes
if ~isempty(p.norm)
    if isequal(p.norm,'tip_norm')
        for i=1:size(ModesU,1)
            fact=1 / ModesU(i,end);
            ModesU(i,:) = ModesU(i,:)*fact;
            ModesV(i,:) = ModesV(i,:)*fact;
            ModesK(i,:) = ModesK(i,:)*fact;
        end
    else
        error('Norm not implemented or incorrect: `%s`',p.norm)
    end
end


