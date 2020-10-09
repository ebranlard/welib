function [T,Y] = fodeRK4(ode,vt,y0,options,acc0)
% Perform Runge-Kutta 4 time integration
%
% Same interface as ode functions of matlab but using time steps prescribed in vt, no adaptative time step!
%
% Author: E. Branlard

%% TEST FUNCTION
if nargin==0
    close all;
    %% Analytical equation
    fPos=@(t)    [ cos(t); sin(t)];
    fVel=@(t)    [-sin(t); cos(t)];
    fAcc=@(t)    [-cos(t);-sin(t)];
    fFirstOrder=@(t,y)[ fVel(t); fAcc(t)];
    vt = 0:0.1:30;
    acc0=fAcc(0);
    [~,Y] = fodeRK4(fFirstOrder,vt,[fPos(0); fVel(0)],[],acc0);
    pos_th=fPos(vt);
    figure,hold all;
    for i=1:length(fPos(0))
        plot(vt,Y(:,i));
        plot(vt,pos_th(i,:),'k.');
        ylim([-1 1]);
    end

    %% Non Linear equation - Parametric excitation
    beta=0.05; sigma=0.01; q=0.2; omega0=1; gamma=-1/6;
    Omega=(2+sigma)*omega0;
    C1=(0.5*q*(Omega/omega0)^2)^2 - (2*beta)^2; % Thomsen 3.140
    a2= sqrt(4/(3*gamma)*(sigma-sqrt(C1) )); % Thomsen 3.139
    vt = 0:0.1:100;
    %fFirstOrder=@(t,y)[y(2) ; -2*beta*omega0*y(2)-(omega0^2-q*Omega^2 *cos(Omega*t))*sin(y(1))];
    fFirstOrder=@(t,y)[y(2) ; -2*beta*omega0*y(2)-(omega0^2-q*Omega^2 *cos(Omega*t))*y(1)-gamma*y(1)^3]; % Thomsen 3.107
    Y0=[pi/6;0];
    tic(); [~,Y]  = fodeRK4(fFirstOrder,vt,Y0); toc();
    tic(); [~,Ym] = ode45  (fFirstOrder,vt,Y0); toc();
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
    vt = 0:0.1:150;
    fFirstOrder=@(t,y)[y(2) ; -2*beta*omega0*y(2)-omega0^2*y(1)-gamma*y(1)^3 + q*cos(Omega*t)]; % Thomsen 3.149
    Y0=[pi/6;0];
    tic(); [~,Y]  = fodeRK4(fFirstOrder,vt,Y0); toc();
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
    fM=@(x) m*eye(2);
    feps=@(x) (x(1)^2+x(2)^2-l0^2)/(2*l0^2);
    fN=@(x) EA*feps(x);
    fg=@(x) [x(1);x(2)]*fN(x)/l0;
    ff=@(x) [m*g;0];
    fAcc=@(x,xp)  fM(x) \ ( ff(x) - fg(x)   ) ;
    fFirstOrder=@(t,y)[y(3:4); fAcc(y(1:2),y(3:4)) ]; 
    %vt = 0:0.0018:6; % NOTE: should be 0.2
    vt = 0:0.002:6; 
    Y0=[0;1.1;0;0];
    tic(); [~,Y]  = fodeRK4(fFirstOrder,vt,Y0); toc();
    opts=odeset('AbsTol',1e-15,'RelTol',1e-7);
    tic(); [~,Ym] = ode45 (fFirstOrder,vt,Y0,opts); toc();
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
    fprintf('DE = %.1e%% (ode45)\n'       ,max(abs(Em/(m*g*l0)-1.65375))/1.65375*100);
    fprintf('DE = %.1e%% (odeRK-cst-dt)\n',max(abs(E /(m*g*l0)-1.65375))/1.65375*100);
    fprintf('Dx = %.1e%% \n'       ,mean(abs((Ym(:,1)-Y(:,1))))/l0*100);
    figure,hold all; title('Nonlinear pendulum');
    plot(vt,Ym(:,1),'k-');
    plot(vt,Ym(:,2),'k-');
    plot(vt,Y(:,1));
    plot(vt,Y(:,2));
%     figure,hold all; title('Nonlinear pendulum - Length and Angle');
%     plot(vt,lm(:,1)/l0,'b-');
%     plot(vt,2*phim/pi,'m--');
    return
end


y0   = y0(:)                 ;
nDOF = length(y0)/2          ;
Y    = zeros(length(vt),nDOF*2);
T    = vt                    ;

pos = y0(1:nDOF);
vel = y0(nDOF+1:end);
if exist('acc0','var') 
    acc = acc0;
else
%     acc = zeros(nDOF,1);
    A0=ode(vt(1), [pos; vel]);
    acc=A0(nDOF+1:end);
end

% --------------------------------------------------------------------------------}
%% --- Default Arguments
% --------------------------------------------------------------------------------{
bProgress=false;
if ~exist('options','var');   options=struct();  end
if isfield(options,'OutputFcn')
    bProgress=true;
    fProgress=options.OutputFcn;
end




% Storing intial value
Y(1,:)  = [pos(:)' vel(:)'];

if bProgress
    fProgress(vt(end),[],'init');
end

for it=2:length(vt)
    t  = vt(it-1)       ;
    dt = vt(it)-vt(it-1);
    %--- Solving for acceleration 
    % 4-step Runge-Kutta-Nystrom time integration of 2nd order differential system integration
    % --- Estimates at t+dt/2
    A=dt/2*acc;
    % Speed changed to v+A, position changed to x+b
    b=dt/2*(vel+A/2);
    B=dt/2 * ode(t+dt/2, [pos+b; vel+A]); 
    B=B(nDOF+1:end);
    % speed changed to V+B
    C=dt/2 * ode(t+dt/2, [pos+b ; vel+B]);
    C=C(nDOF+1:end);
    % --- Estimates at t+dt
    % speed changed to v+2*C, position changed to x+d
    d=dt*(vel+C);
    D=dt/2 * ode(t+dt , [pos+d ; vel+2*C]);
    D=D(nDOF+1:end);
    % final estimate
    pos = pos+dt*(vel+1/3*(A+B+C)) ;
    vel = vel+1/3*(A+2*B+2*C+D)  ;
    acc = ode(t+dt , [pos ; vel]);
    acc = acc(nDOF+1:end);
    % Storing 
    Y(it,:)  = [pos(:)' vel(:)'];
    if bProgress
        fProgress(t+dt,[],[]);
    end
end
if bProgress
    fProgress(vt(end),[],'done');
end
