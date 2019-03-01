function [ P PS SS Cp u v xg yg chord X_profile Y_profile] = fProfileVanDeVooren( t_rel,tau,chord,n,flagMethod,varargin )
% returns Van De Vooren profiles coordinates and optionally some aero data 
% --- Input ------------------
% tau     = 10   ;  % Tail Angle in deg !!!
% t_rel = 0.15 ;  % approx t/c (approx ymax_SS-ymin_PS), later thethickness parameter e is defined as t_rel/2
% alpha = 0    ; 
% n       = 201  ;  % Number of points along mapped foil surface including twice the TE
% out_pts = 36   ;  % Total number of X-Y output points (ODD)
% --- Output -------------
% P=[x y]: normalized coordinates of the airfoil
% PS idem but pressure side
% SS idem but suction side
% Cp : pressure coeff
% u,v,xg,yg : velocity field around profile and grid points (polar grid>use TriScatteredInterp to go back to a cartesian grid)
% chord
%% Input checks
if nargin>6
    bComputeAero=1;
    U0=varargin{1};
    alpha=varargin{2}; % !!!!!!!!!!!!!!!!!!!!!! alpha in degrees!!!
    %     if nargin>6
    %         % velocity computation on polar Grid
    %         bComputeGridVelocity=1;
    %         vr=varargin{7}; % radial vector 
    %         vt=varargin{8}; % tangential vector
    %     else
    %         bComputeGridVelocity=0;
    %     end
else
    bComputeAero=0;
end
% Optional outputs
Cp=[]; u=[]; v=[]; xg=[]; yg=[]; 

P=[];
PS=[];
SS=[];


%% Calculated Properties
thickness=t_rel/2;
lambda = 2-tau/180; % lambda= 2-tau_rad/pi
if chord~=1
    error('I have to rethink stuff again, for now chord==1')
end

l = chord/2; % not really needed, but for consistency
a = chord*(1+thickness)^(lambda-1)*2^(-lambda); 



% ---------------------------------------------------------------------------
%% Method1 : method involving the function fConformalMap and fUi_Cylinder2D
% ---------------------------------------------------------------------------
if flagMethod==1
    % k

    if ~exist('fConformalVanDeVooren','file')
        require('POTFLOW','v00')
    end

    % Circle
    vtheta_circ = 0:2*pi/n_pts:2*pi-pi/n_pts; %Defines theta incremented 0->2*pi
    z_circ= a*exp(i*vtheta_circ);
    % Karman Trefftz Conformal map
    [ Z_profile dZdz] = fConformalMapVanDeVooren( z_circ, a, lambda,thickness,l  );
    X_profile=real(Z_profile);
    Y_profile=imag(Z_profile);

    %% Aero computation if required -> this is just minimum if needed, see dedicated script for more..
    if bComputeAero
        Gamma = 4*pi*a*U0*sind(-alpha); %Kutta Condition (beta=0)
        % Velocity at circle
        [ u_circ, v_circ ] = fUi_Cylinder2D(real(z_circ),imag(z_circ), xc,yc,rc, U0, alpha, Gamma  , false);
        % Velocities -Cp on surface
        W_circ  = (u_circ-i*v_circ)./dZdz    ;  % [u-iv]_Z = [u-iv]_z/DZ/Dz
        U_circ  = real(W_circ)           ;  % X velocity in Zeta-plane
        V_circ  = -imag(W_circ)          ;  % Y velocity in Zeta-plane
        Q  = sqrt(U_circ.^2 + V_circ.^2) ;  % Z-plane Velocity Magnitude
        Cp = 1-(Q./U0).^2      ;  % Z-plane pressure coefficient

        %         %% direct transform,polar coord grid
        %         if bComputeGridVelocity
        %             [r,theta] = meshgrid(vr+rc,vt); %Create location mesh
        %             xg=r.*cos(theta)+real(z0);
        %             yg=r.*sin(theta)+imag(z0);
        %             [ u, v ]  = fUi_Cylinder2D(xg,yg, xc,yc,rc, U0, alpha, Gamma  ,false);
        %         end
    end


    % ---------------------------------------------------------------------------
    %% Method2 : standalone method inspired by Katz Prog1
    % ---------------------------------------------------------------------------
elseif flagMethod==2
    % based on Katz program 1 
    e=thickness;
    ak=lambda;
    a = chord*(1+thickness)^(lambda-1)*2^(-lambda); 
    l = chord/2; % not really needed, but for consistency

    theta=linspace(0,2*pi,n);

    r1=sqrt((a*(cos(theta)-1)).^2+(a*sin(theta)).^2) ;
    r2=sqrt((a*(cos(theta)-e)).^2+(a*sin(theta)).^2) ;

    % default: % (simplified by a)
    theta1=atan(sin(theta)./(cos(theta)-1))+pi;
    theta2=atan(sin(theta)./((cos(theta)-e)));

    % exceptions and co
    theta1(theta==0)=pi/2;
    I=cos(theta)-e<0 & sin(theta)>0;
    theta2(I)=theta2(I)+pi;
    I=(cos(theta)-e<0 & sin(theta)<0);
    theta2(I)=theta2(I)+pi;
    I=(cos(theta)-e>0 & sin(theta)<0);
    theta2(I)=theta2(I)+2*pi;

    com1=( (r1.^ak)./(r2.^(ak-1)) )./((cos((ak-1) *theta2)).^2+(sin((ak-1)*theta2)).^2);
    X_profile=com1.*(cos(ak*theta1).*cos((ak-1)*theta2) +sin(ak*theta1).*sin((ak-1)*theta2))+l;
    Y_profile=com1.*(sin(ak*theta1).*cos((ak-1)*theta2) -cos(ak*theta1).*sin((ak-1)*theta2));


    % C1=(cos(ak*theta2)).^2+(sin(ak*theta2)).^2

    if bComputeAero
        alphar=alpha*pi/180;

        A1=cos((ak-1)*theta1).*cos(ak*theta2)+sin((ak-1)*theta1).*sin(ak*theta2);
        B1=sin((ak-1)*theta1).*cos(ak*theta2)-cos((ak-1)*theta1).*sin(ak*theta2);

        D0=a*(1-ak+ak*e);
        D1=A1.*(a*cos(theta)-D0)-B1*a.*sin(theta) ;
        D2=A1*a.*sin(theta)+B1.*(a*cos(theta)-D0) ;
        
        com2=2*(r2.^ak)./(r1.^(ak-1)) .* (sin(alphar)-sin(alphar-theta))./(D1.^2+D2.^2) ;
        vx=  D1.*sin(theta)+D2.*cos(theta) ;
        vy=-(D1.*cos(theta)-D2.*sin(theta)) ;
        Cp=1-com2.^2.*(vx.^2+vy.^2);
    end
end


[P PS SS TE chord]=fProfileStandardize(X_profile,Y_profile,1);

