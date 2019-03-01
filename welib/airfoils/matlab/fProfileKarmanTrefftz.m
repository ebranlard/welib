function [ P PS SS Cp u v xg yg chord X_profile Y_profile IPin IPout ] = fProfileKarmanTrefftz( xc,yc,tau,n,varargin )
% returns Van Karman profiles and optionally some aero data 
% --- Input ------------------
% xc      = -0.2 ;  % Circle Center Location (<0)
% yc      = -0.1 ;  % Circle Center (>0 adds + camber)
% tau     = 10   ;  % Tail Angle in deg !!!
% n       = 201  ;  % Number of points along mapped foil surface
% out_pts = 36   ;  % Total number of X-Y output points (ODD)
% --- Output -------------
% P=[x y]: normalized coordinates of the airfoil
% PS idem but pressure side
% SS idem but suction side
% Cp : pressure coeff
% u,v,xg,yg : velocity field around profile and grid points (polar grid>use TriScatteredInterp to go back to a cartesian grid)
% chord
if ~exist('fConformalMapKarmanTrefftz','file')
    require('POTFLOW','v00')
end


if nargin>4
    bComputeAero=1;
    U0=varargin{1};
    alpha=varargin{2}; % deg!!!
    if nargin>6
        % velocity computation on polar Grid
        bComputeGridVelocity=1;
        vr=varargin{3}; % radial vector 
        vt=varargin{4}; % tangential vector
    else
        bComputeGridVelocity=0;
    end
else
    bComputeAero=0;
end
% Optional outputs
Cp=[]; u=[]; v=[]; xg=[]; yg=[]; 

% Param and Properties
a      = 1.0                   ;  % x-intersect
rc     = sqrt((a-xc)^2 + yc^2) ;  % radius of circle
beta   = asin(-yc/(rc))        ;  % Angle to rear stagnation point
lambda = 2-tau/180             ; 

% Circle
vtheta_circ = 0:2*pi/n:2*pi-pi/n; %Defines theta incremented 0->2*pi
z0=(xc + i*yc); % center of circle
z_circ= z0 + rc*exp(i*vtheta_circ);

% Karman-Trefftz Conformal map
[ Z_profile dZdz] = fConformalMapKarmanTrefftz( z_circ, a, lambda  );
X_profile=real(Z_profile(:));
Y_profile=imag(Z_profile(:));


% figure, 
% fPlotDegrad(xps,yps,2);
% fPlotDegrad(xss,yss,1);

%%


%% Aero computation if required -> this is just minimum if needed, see dedicated script for more..
if bComputeAero
    Gamma = 4*pi*rc*U0*sin(beta-alpha*pi/180) % from Kutta condition 
    % Velocity at circle
    [ u_circ, v_circ ] = fUi_Cylinder2D(real(z_circ),imag(z_circ), xc,yc,rc, U0, alpha, Gamma  , false);
    % Velocities -Cp on surface
    W_circ  = (u_circ-i*v_circ)./dZdz    ;  % [u-iv]_Z = [u-iv]_z/DZ/Dz
    U_circ  = real(W_circ)           ;  % X velocity in Zeta-plane
    V_circ  = -imag(W_circ)          ;  % Y velocity in Zeta-plane
    Q  = sqrt(U_circ.^2 + V_circ.^2) ;  % Z-plane Velocity Magnitude
    Cp = 1-(Q./U0).^2      ;  % Z-plane pressure coefficient

%     xCP=(X_profile-min(X_profile))./(max(X_profile)-min(X_profile));


    %% direct transform,polar coord grid
    if bComputeGridVelocity
        [r,theta] = meshgrid(vr+rc,vt); %Create location mesh
        xg=r.*cos(theta)+real(z0);
        yg=r.*sin(theta)+imag(z0);
        [ u, v ]  = fUi_Cylinder2D(xg,yg, xc,yc,rc, U0, alpha, Gamma  ,false);
    end
end


%% Standardizing coordinates
[P PS SS TE chord IPin IPout]=fProfileStandardize(X_profile,Y_profile,[]);



