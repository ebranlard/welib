function [T,Y,xpp] = fodeNewmarkNonLinear(fMDKgf,vt,y0,options)
%
% Solving the equation:
%   xddot = M^-1 (g) (=) M^{-1}( F - D xdot - Kx ) 
%
% INPUTS:
%    fMDKgf=@function [M,D,K,g,f]=fMDKGf(t,x,xp) 
%
% Similar interface as ode functions of matlab
%
% Author: E. Branlard
%
%% TEST FUNCTION
if nargin==0
    close all; clear all;

    %% Pendulum with elastic bar - Krenk
    m=1; l0=1; g=10; EA=3000;  
    fM   = @(x) m*eye(2);
    fD   = @(x) zeros(2,2);
    feps = @(x) (x(1)^2+x(2)^2-l0^2)/(2*l0^2);
    fN   = @(x) EA*feps(x);
    fg   = @(x,xp) [x(1) ; x(2)]*fN(x)/l0;
    ff   = @(x)    [m*g  ; 0]            ;
    fK   = @(x) eye(2)* fN(x)/l0 + EA/l0^3*x*x';
    fAcc = @(x,xp)  fM(x) \ ( ff(x) - fg(x)   );
    options.gamma=1/2;
    options.beta=1/4;
    options.epsilon_u=1e-6*m*g;
    options.epsilon_r=1e-6*l0;
    fFirstOrder=@(t,y)[y(3:4); fAcc(y(1:2),y(3:4)) ]; 
    fMDKgf=@(t,x,xp) deal( fM(x), fD(x), fK(x), fg(x,xp), ff(x)) ;
    vt = 0:0.02:6; % NOTE: should be 0.02
    Y0=[0;1.1;0;0];
    tic(); [~,Y]  = fodeNewmarkNonLinear(fMDKgf,vt,Y0,options); toc();
    options=odeset('AbsTol',1e-15,'RelTol',1e-7);
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0,options); toc();
    lm   = sqrt(Ym(:,1).^2+Ym(:,2).^2)          ;
    epsm = (lm.^2-l0^2)/(2*l0^2);
    Um   = 0.5*l0*EA*epsm.^2 - m*g*Ym(:,1)      ;
    Tm   = 0.5*m*(Ym(:,3).^2+Ym(:,4).^2)        ;
    Em=Um+Tm;
    l   = sqrt(Y(:,1).^2+Y(:,2).^2)          ;
    eps = (l.^2-l0^2)/(2*l0^2);
    U   = 0.5*l0*EA*eps.^2 - m*g*Y(:,1)      ;
    T   = 0.5*m*(Y(:,3).^2+Y(:,4).^2)        ;
    E=U+T;
    figure,hold all; title('Nonlinear pendulum Total Energy');
    plot(vt,Em/(m*g*l0),'k-');
    plot(vt,E /(m*g*l0),'b-');
    fprintf('DE = %.1e%% (ode45)\n'     ,max(abs(Em/(m*g*l0)-1.65375))/1.65375*100);
    fprintf('DE = %.1e%% (NL-Newmark)\n',max(abs(E /(m*g*l0)-1.65375))/1.65375*100);
    fprintf('Dx = %.1e%% \n'            ,mean(abs((Ym(:,1)-Y(:,1))))/l0*100);
    figure,hold all; title('Nonlinear pendulum');
    plot(vt,Ym(:,1),'k-');
    plot(vt,Ym(:,2),'k-');
    plot(vt,Y(:,1));
    plot(vt,Y(:,2));

    %% Pendulum with elastic bar - Krenk - No K C
    m=1; l0=1; g=10; EA=3000;  
    fM   = @(x) m*eye(2);
    fD   = @(x) zeros(2,2);
    fK   = @(x) zeros(2,2);
    feps = @(x) (x(1)^2+x(2)^2-l0^2)/(2*l0^2);
    fN   = @(x) EA*feps(x);
    fg   = @(x,xp) [x(1) ; x(2)]*fN(x)/l0;
    ff   = @(x)    [m*g  ; 0]            ;
    fAcc = @(x,xp)  fM(x) \ ( ff(x) - fg(x)   );
    options.gamma=1/2;
    options.beta=1/4;
    options.epsilon_u=1e-6*m*g;
    options.epsilon_r=1e-6*l0;
    fFirstOrder=@(t,y)[y(3:4); fAcc(y(1:2),y(3:4)) ]; 
    fMDKgf=@(t,x,xp) deal( fM(x), fD(x), fK(x), fg(x,xp), ff(x)) ;
    vt = 0:0.02:6; % NOTE: should be 0.02
    Y0=[0;1.1;0;0];
    tic(); [~,Y]  = fodeNewmarkNonLinear(fMDKgf,vt,Y0,options); toc();
    options=odeset('AbsTol',1e-15,'RelTol',1e-7);
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0,options); toc();
    lm   = sqrt(Ym(:,1).^2+Ym(:,2).^2)          ;
    epsm = (lm.^2-l0^2)/(2*l0^2);
    Um   = 0.5*l0*EA*epsm.^2 - m*g*Ym(:,1)      ;
    Tm   = 0.5*m*(Ym(:,3).^2+Ym(:,4).^2)        ;
    Em=Um+Tm;
    l   = sqrt(Y(:,1).^2+Y(:,2).^2)          ;
    eps = (l.^2-l0^2)/(2*l0^2);
    U   = 0.5*l0*EA*eps.^2 - m*g*Y(:,1)      ;
    T   = 0.5*m*(Y(:,3).^2+Y(:,4).^2)        ;
    E=U+T;
    figure,hold all; title('Nonlinear pendulum Total Energy');
    plot(vt,Em/(m*g*l0),'k-');
    plot(vt,E /(m*g*l0),'b-');
    ylim([1.4 2])
    fprintf('DE = %.1e%% (ode45)\n'     ,max(abs(Em/(m*g*l0)-1.65375))/1.65375*100);
    fprintf('DE = %.1e%% (NL-Newmark)\n',max(abs(E /(m*g*l0)-1.65375))/1.65375*100);
    fprintf('Dx = %.1e%% \n'            ,mean(abs((Ym(:,1)-Y(:,1))))/l0*100);
    figure,hold all; title('Nonlinear pendulum');
    plot(vt,Ym(:,1),'k-');
    plot(vt,Ym(:,2),'k-');
    plot(vt,Y(:,1));
    plot(vt,Y(:,2));


    return
end


%% Init of output variables
y0   = y0(:)                 ;
nDOF = length(y0)/2          ;
x    = zeros(nDOF,length(vt));
xp   = zeros(nDOF,length(vt));
xpp  = zeros(nDOF,length(vt));
T    = vt                    ;

% --------------------------------------------------------------------------------}
%% --- Default Arguments
% --------------------------------------------------------------------------------{
bProgress=false;
if ~exist('options','var');   options=struct();  end
if ~isfield(options,'gamma'); options.gamma=1/2; end
if ~isfield(options,'beta');  options.beta=1/4;  end
if ~isfield(options,'epsilon_u');  options.epsilon_u=1e-3;  end
if ~isfield(options,'epsilon_r');  options.epsilon_r=1e-3;  end
if isfield(options,'OutputFcn')
    if ~isempty(options.OutputFcn)
        bProgress=true;
        fProgress=options.OutputFcn;
    end
end

nItMax=20;

gamma = options.gamma;
beta  = options.beta ;


% Initial value - Solving done to find initial acceleration
x (:,1)= y0(1:nDOF    );
xp(:,1)= y0(nDOF+1:end);
[M,~,~,g,f] = fMDKgf(vt(1), x(:,1), xp(:,1)  );
xpp(:,1) = M \ (f-g) ;

if bProgress
    fProgress(vt(end),[],'init');
end


%% Time loop
nItTot=0;
for it=1:length(vt)-1
    t  = vt(it+1)       ;
    dt = vt(it+1)-vt(it);

    % ---  Krenk p.305
    % Prediction step
    xpps=  xpp(:,it);
    xps =  xp (:,it)  + dt*xpp(:,it);
    xs  =  x  (:,it)  + dt*xp (:,it) + 0.5*dt^2*xpp(:,it);
    % System evaluation
    [M,C,K,g,f] = fMDKgf(t, xs, xps);
    % Residual calculation
    r = f - M*xpps - g;
    dx=10*options.epsilon_u; % to avoid success criteria 
    nIt=0;
    while (norm(r)>options.epsilon_r && norm(dx)>options.epsilon_u && nIt<nItMax)
        nIt=nIt+1;
        % Correction step 
        Ks = K + gamma/(beta*dt)*C + 1/(beta*dt^2)*M;
        dx = Ks\r;
        xs   = xs   +                 dx;
        xps  = xps  + gamma/(beta*dt)*dx;
        xpps = xpps + 1/(beta*dt^2)  *dx;
        % System evaluation
        [~,~,~,g,f] = fMDKgf(t, xs, xps);
        % Residual calculation
        r = f - M*xpps - g;
    end
    nItTot=nItTot+nIt;
    if nIt>=nItMax
        fprintf('Max iteration reached at t=%.3f\n',t);
    end
    % Update
    x  (:,it+1) = xs  ;
    xp (:,it+1) = xps ;
    xpp(:,it+1) = xpps;
    if bProgress
        fProgress(t,[],[]);
    end
end
if bProgress
    fProgress(vt(end),[],'done');
end
fprintf('Average number of iterations: %.1f\n',nItTot/(length(vt)-1));
%% "Matlab ode format"
Y=[x' xp'];
