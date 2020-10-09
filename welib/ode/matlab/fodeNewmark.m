function [T,Y,xpp] = fodeNewmark(fMDKR,vt,y0,options)
% Perform Newmark time integration
%
% Similar interface as ode functions of matlab but using time steps prescribed in vt, no adaptative time step!
% Also the function handle fMDKR represents a function that returns M,K,D,R
%
% Author: E. Branlard
%
%% TEST FUNCTION
if nargin==0
    close all;

    %% Forced harmonic vibration  - No Phase
    k = 5.7e+06; m = 2.2e+05; zeta = 0.5/(2*pi); omega0 = sqrt(k/m); Omega=0.99*omega0; F0= 1e6;
    fMDKR=@(t,x,xp) deal(m,2*m*omega0*zeta,k,F0*sin(Omega*t));
    vt = 0:0.10:10;
    tic(); [~,Y] = fodeNewmark(fMDKR,vt,[0;0]); toc();
    frat =  Omega/omega0;
    A_th   = (F0/k) * 1./sqrt( (1-frat.^2).^2 + (2*zeta*frat).^2 )                                            ;
    phi_th = atan(2*zeta*frat/(1-frat^2))                                                                     ;
    psi    = atan( sqrt(1-zeta^2) / (zeta-Omega/omega0*cos(phi_th)/sin(phi_th)))                              ;
    A      = A_th* sin(phi_th)/sin(psi)                                                                       ;
    x_th   = A_th * sin(Omega*vt - phi_th) + A*exp(-zeta*omega0*vt) .* sin(omega0 * sqrt(1-zeta^2) * vt + psi);

    figure,hold all;
    plot(vt,Y(:,1));
    plot(vt,x_th,'k.');

    %% Forced harmonic vibration  - With phase and Y0
    k = 5.7e+06; m = 2.2e+05; zeta = 0.5/(2*pi); omega0 = sqrt(k/m); Omega=0.99*omega0; F0= 1e6; Phi= pi/2;
    fFirstOrder=@(t,x)[x(2) ; -2*omega0*zeta*x(2)-k/m*x(1)+F0/m*sin(Omega*t+Phi) ]; 
    fMDKR=@(t,x,xp) deal(m,2*m*omega0*zeta,k,F0*sin(Omega*t+Phi));
    vt = 0:0.10:10;
    Y0 = [pi/6; pi/10];
    tic(); [~,Y]  = fodeNewmark(fMDKR,vt,Y0)    ; toc();
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0); toc();
    figure,hold all;
    plot(vt,Y(:,1));
    plot(vt,Ym(:,1),'k.');

    %% Non Linear equation - Parametric excitation
    beta=0.05; sigma=0.01; q=0.2; omega0=1; gamma=-1/6;
    Omega=(2+sigma)*omega0;
    C1=(0.5*q*(Omega/omega0)^2)^2 - (2*beta)^2; % Thomsen 3.140
    a2= sqrt(4/(3*gamma)*(sigma-sqrt(C1) )); % Thomsen 3.139
    vt = 0:0.20:100;
    fFirstOrder=@(t,y)[y(2) ; -2*beta*omega0*y(2)-(omega0^2-q*Omega^2 *cos(Omega*t))*y(1)-gamma*y(1)^3]; % Thomsen 3.107
    fMDKR=@(t,x,xp)deal(1, 2*beta*omega0, omega0^2-q*Omega^2 *cos(Omega*t),-gamma*x^3); % Thomsen 3.107    
%     fMDKR=@(t,x,xp)deal(1, 2*beta*omega0, omega0^2,-gamma*x^3+q*Omega^2 *cos(Omega*t)*x); % Thomsen 3.107    
    Y0=[pi/6; 0];
    tic(); [~,Y]  = fodeNewmark(fMDKR,vt,Y0)    ; toc();
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0); toc();
    figure,hold all;
    plot(vt,Y(:,1));
    plot(vt,Ym(:,1),'k.');
    plot(vt,vt*0+a2,'k--');
    plot(vt,vt*0-a2,'k--');
    ylim([-a2 a2]*1.5);


    %% Non Linear equation - Duffing equation
    beta=0.05; a=0.9; q=0.2; omega0=1; gamma=0.5;
    %OmegaM=omega0*(1+3*gamma/(8*omega0^2)*a^2 - sqrt( (q/(2*omega0^2*a))^2-beta^2 )); % Thomsen 3.194
    OmegaP=omega0*(1+3*gamma/(8*omega0^2)*a^2 + sqrt( (q/(2*omega0^2*a))^2-beta^2 )); % Thomsen 3.194
    Omega=OmegaP; 
    vt = 0:0.2:150; % NOTE: with a smaller time step is works... But... 
    fFirstOrder=@(t,y)[y(2) ; -2*beta*omega0*y(2)-omega0^2*y(1)-gamma*y(1)^3 + q*cos(Omega*t)]; % Thomsen 3.149
    fMDKR=@(t,x,xp)deal(1, 2*beta*omega0, omega0^2,-gamma*x^3+q*cos(Omega*t)); % Thomsen 3.149    
    Y0=[pi/6;0];
    tic(); [~,Y]  = fodeNewmark(fMDKR,vt,Y0)    ; toc();
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0); toc();
    figure,hold all; title('Duffing equation');
    plot(vt,Y(:,1));
    plot(vt,Ym(:,1),'k.');
    aa=a+gamma/(32*omega0^2)*a^3;
    plot(vt,vt*0+aa,'k--');
    plot(vt,vt*0-aa,'k--');
    ylim([-aa aa]*1.5);


    %% Pendulum with elastic bar - Krenk 
    m=1; l0=1; g=10; EA=3000;  
    fM   = @(x) m*eye(2);
    fD   = @(x) zeros(2,2);
    fK   = @(x) zeros(2,2);
    feps = @(x) (x(1)^2+x(2)^2-l0^2)/(2*l0^2);
    fN   = @(x) EA*feps(x);
    fg   = @(x,xp) [x(1) ; x(2)]*fN(x)/l0;
    ff   = @(x)    [m*g  ; 0]            ;
    fAcc = @(x,xp)  fM(x) \ ( ff(x) - fg(x)   );
    fFirstOrder=@(t,y)[y(3:4); fAcc(y(1:2),y(3:4)) ]; 
    fMDKR=@(t,x,xp) deal( fM(x), fD(x), fK(x), ff(x)-fg(x,xp)) ;
    vt = 0:0.02:6; % NOTE: should be 0.02
    Y0=[0;1.1;0;0];
    tic(); [~,Y]  = fodeNewmark(fMDKR,vt,Y0); toc();
    opts=odeset('AbsTol',1e-15,'RelTol',1e-7);
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0,opts); toc();
    lm   = sqrt(Ym(:,1).^2+Ym(:,2).^2)          ;
    %phim = atan(Ym(:,2)./Ym(:,1))               ;
    epsm = (lm.^2-l0^2)/(2*l0^2);
    Um   = 0.5*l0*EA*epsm.^2 - m*g*Ym(:,1)      ;
    Tm   = 0.5*m*(Ym(:,3).^2+Ym(:,4).^2)        ;
    Em=Um+Tm;
    l   = sqrt(Y(:,1).^2+Y(:,2).^2)          ;
    %phi = atan(Y(:,2)./Y(:,1))               ;
    eps = (l.^2-l0^2)/(2*l0^2);
    U   = 0.5*l0*EA*eps.^2 - m*g*Y(:,1)      ;
    T   = 0.5*m*(Y(:,3).^2+Y(:,4).^2)        ;
    E=U+T;
    figure,hold all; title('Nonlinear pendulum Total Energy');
    plot(vt,Em/(m*g*l0),'k-');
    plot(vt,E /(m*g*l0),'b-');
    fprintf('DE = %.1e%% (ode45)\n'  ,max(abs(Em/(m*g*l0)-1.65375))/1.65375*100);
    fprintf('DE = %.1e%% (Newmark)\n',max(abs(E /(m*g*l0)-1.65375))/1.65375*100);
    fprintf('Dx = %.1e%% \n'         ,mean(abs((Ym(:,1)-Y(:,1))))/l0*100);
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
if isfield(options,'OutputFcn')
    if ~isempty(options.OutputFcn)
        bProgress=true;
        fProgress=options.OutputFcn;
    end
end
gamma = options.gamma;
beta  = options.beta ;


% Initial value - Solving done to find initial acceleration
x (:,1)= y0(1:nDOF    );
xp(:,1)= y0(nDOF+1:end);
[M,C,K,R0] = fMDKR(vt(1), x(:,1), xp(:,1)  );
xpp(:,1) = M \ (R0- C*xp(:,1) - K*x(:,1)) ;
if bProgress
    fProgress(vt(end),[],'init');
end



%% Time loop
for it=1:length(vt)-1
    t  = vt(it+1)       ;
    dt = vt(it+1)-vt(it);

    % --- Formulation 1 - Krenk p.295
    % Prediction step
    xps =  xp (:,it)                 +  (1-gamma)*dt  *xpp(:,it);
    xs  =  x  (:,it)  + dt* xp(:,it) + (1/2-beta)*dt^2*xpp(:,it);
    % Correction step and update
    [M,C,K,Rts] = fMDKR(t, xs, xps  );
    Ms = M + gamma*dt*C + beta* dt^2 * K;
    xpp(:,it+1) = Ms\(Rts- C*xps - K*xs);
    xp (:,it+1) = xps + gamma*dt  *xpp(:,it+1);
    x  (:,it+1) = xs  + beta *dt^2*xpp(:,it+1);

    % --- Formulation 2
%     % Getting matrices and force
%     [M,C,K,Rt] = fMDKR(t, x(:,it), xp(:,it)  );
% 
%     % Constants used in Newmark's integration
%     a1 = gamma/(beta*dt);
%     a2 = 1/(beta*dt^2);
%     a3 = 1/(beta*dt);
%     a4 = gamma/beta;
%     a5 = 1/(2*beta);
%     a6 = (gamma/(2*beta)-1)*dt;
% 
%     A   = K + a1*C + a2*M    ; % Needs to be inverted to find dx
%     RHS = (Rt-R0) + (a3*M + a4*C)*xp(:,it) + (a5*M + a6*C)*xpp(:,it);
%     dx  = A \ RHS                                  ; % Getting dx
% 
%     x  (:,it+1) = x  (:,it) +    dx;
%     xp (:,it+1) = xp (:,it) + a1*dx - a4*xp(:,it)- a6*xpp(:,it);
%     xpp(:,it+1) = xpp(:,it) + a2*dx - a3*xp(:,it)- a5*xpp(:,it);
%     R0=Rt;
    if bProgress
        fProgress(t,[],[]);
    end
end
if bProgress
    fProgress(vt(end),[],'done');
end

%% "Matlab ode format"
Y=[x' xp'];
