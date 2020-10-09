function [T,Y] = fodeKrenk(fMDKgKcGF,vt,y0,opts)
% Implement Krenk numerial integration scheme
%
% [M,D,Kgm,Kcm,gm,fm] = fMDKgKcGF(t, um, vn);
%
% Apart from that, similar interface as ode functions of matlab 
%
% Author: E. Branlard
%
%% TEST FUNCTION
if nargin==0
    close all; clear all;

    %% Pendulum with elastic bar - Krenk 
    % Note if stiffness is 0 then it's more or less like Newmark-NL
    m=1; l0=1; g=10; EA=3000;  
    fM   = @(x) m*eye(2);
    fD   = @(x) zeros(2,2);
    feps = @(x) (x(1)^2+x(2)^2-l0^2)/(2*l0^2);
    fN   = @(x) EA*feps(x);
    fg   = @(x,xp) [x(1) ; x(2)]*fN(x)/l0;
    ff   = @(x)    [m*g  ; 0]            ;
    fKg  = @(x) eye(2)* fN(x)/l0;
    fKc  = @(x) EA/l0^3*x*x';
    fAcc = @(x,xp)  fM(x) \ ( ff(x) - fg(x)   );
    opts.epsilon_u=1e-8*m*g;
    opts.epsilon_r=1e-8*l0;
    opts.epsilon_g=1e-6*l0; %TODO
    fFirstOrder=@(t,y)[y(3:4); fAcc(y(1:2),y(3:4)) ]; 
    fMDKgKcGF=@(t,x,xp) deal( fM(x), fD(x), fKg(x), fKc(x), fg(x,xp), ff(x)) ;
    vt = 0:0.002:6; 
    Y0=[0;1.1;0;0];
    tic(); [~,Y]  = fodeKrenk(fMDKgKcGF,vt,Y0,opts); toc();
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
    fprintf('DE = %.1e%% (ode45)\n',max(abs(Em/(m*g*l0)-1.65375))/1.65375*100);
    fprintf('DE = %.1e%% (Krenk)\n',max(abs(E /(m*g*l0)-1.65375))/1.65375*100);
    fprintf('Dx = %.1e%% \n'       ,mean(abs((Ym(:,1)-Y(:,1))))/l0*100);
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
T    = vt                    ;

% --------------------------------------------------------------------------------}
%% --- Default Arguments
% --------------------------------------------------------------------------------{
if ~exist('opts','var') 
    opts.epsilon_u=1e-3;
    opts.epsilon_r=1e-3;
end

nItMax=20;

% Initial value 
x (:,1)= y0(1:nDOF    );
xp(:,1)= y0(nDOF+1:end);
[M,C,Kgn,Kcn,gn,fn] = fMDKgKcGF(vt(1), x(:,1), xp(:,1)  );
vn= xp(:,1);
un= x (:,1);

%% Time loop
nItTot=0;
for it=1:length(vt)-1
    t  = vt(it+1)       ;
    dt = vt(it+1)-vt(it);

    % ---  Krenk p.318
    % Prediction step
    du = dt * vn;
    bNotConverged=true;
    nIt=0;
    while bNotConverged
        % Residual calculation
        um = un+du;
        [~,~,Kgm,Kcm,gm,fm] = fMDKgKcGF(t, um, vn  );
        dKg= Kgm-Kgn;
        r  = fm+fn - (gm+gn) + 4/dt*M*vn  - (  (2/dt)^2*M + (2/dt)*C - 0.5*dKg ) *du;
        % Displacement sub-increment
        Kcs=0.5*(Kcm+Kcn);
        Kgs=0.5*(Kgm+Kgn);
        Ks= Kcs+Kgs+(2/dt)^2*M +(2/dt)*C;
        ddu=Ks\r;
        du=du+ddu;

        nIt=nIt+1;
        bNotConverged= (norm(r)>opts.epsilon_r && norm(ddu)>opts.epsilon_u && nIt<nItMax);
    end
    nItTot=nItTot+nIt;
    if nIt>=nItMax
        fprintf('Max iteration reached at t=%.3f\n',t);
    end
    % Update
    um = un+du;
    vm = (2/dt)*du - vn;
    x  (:,it+1) = um ;
    xp (:,it+1) = vm ;

    % Preparing next step
    un  = um ;
    vn  = vm ;
    fn  = fm ;
    gn  = gm ;
    Kgn = Kgm;
    Kcn = Kcm;

end
fprintf('Average number of iterations: %.1f\n',nItTot/(length(vt)-1));
%% "Matlab ode format"
Y=[x' xp'];
