function [ G ] = fGoldsteinFactor_Matlab( l_bar,B,vx )
    %l_bar=h/(2piR)
    %l_bar=l/R;    
    if max(vx)>1
        error('Vector should be normalized to one - dimensionless');
%         disp('ERROR: Vector should be normalized to one - dimensionless');
    end
    if vx(1)>0.01
%         warning('You are probably thinking about a hub - it used to be doable - check it')
        error('You are probably thinking about a hub - it used to be doable - check it');
    end
%     if(B>100)
% %         warning('Too many blades for Goldstein theory');
% %         G=vx(:)'+1;        
%         error('Too many blades for Goldstein theory');
%     else
    vx=vx(:)';
    %% a not pretty mainupaltion to ensure computation is done within 0 1, since the algorithm below is meant for it (kind of, actually not, there is something about hub radius
    drop0=false;
    drop1=false;
    if(vx(1)~=0)
        drop0=true;
        vx=[0 vx];
    end
    if(vx(end)~=1)
        drop1=true;
        vx=[vx 1];
    end
    R=1; % vx should go from 0 to 1
    w=1; % circulation is linear in w
    G = fGoldsteinCirculation(l_bar*R,w,R,B,vx)*B/(2*pi*l_bar*R*w); 
    if(drop1)
        G=G(1:end-1);
    end
    if(drop0)
        G=G(2:end);
    end
    %
    G(G<0)=0;
    G(isnan(G))=0;
    G=G(:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [ GammaGoldstein] = fGoldsteinCirculation(l_bar,w,R,B,Vr)
drop=0;
if(Vr(1)==0)
   Vr=Vr(2:end);
   drop=1;
end
%% Inititalization
% Position of the helical vortices
N=length(Vr);
% Position of the Control Points
% The control points are placed between the vorticies
% They are used to compute velocities on the vortex sheet
VrCP=(3/(2*N):1/N:1)*R;   % [m] 

%% Calculation of matrices of influence at control points for a helical band
% Circulation is unitary
A_x=zeros(N,N); %
for j=1:length(Vr)  % loop on vorticies
    for i=1:length(VrCP) % loop on Control points
        a=Vr(j);  %[m] radius of helical vortex
        r=VrCP(i);%[m] radius of control point
        A_x(i,j)=fUi_HelixNTheory(1,r,a ,l_bar*R,0,B);            
    end
end

% Condition sum of gamma=0
A_x(end,:)=1;

% Calculation of boundary conditions values on the vortex sheet
% Boundary conditions are defined in equation 27, we apply then on the
% line made by the control points 
U_x0=zeros(N,1);
for i=1:(N-1)
     r=VrCP(i);      %[m]   
     U_x0(i)= w*1/(1+l_bar^2/(r/R)^2);
end
%% Solving for Gamma
Gamma=A_x\U_x0; 
GammaGoldstein=cumsum(Gamma);
if(drop)
    GammaGoldstein=[0;GammaGoldstein];
end
end

function [ ui ] = fUi_HelixNTheory(Gamma,r,r0 ,l,psih,nB);
% this function takes a real l not l_bar
fact=1; % should be one half for lifting line but we apply BC with far wake induced velocities
sign=-1; % for wt, should not matter here since we are on the helix and don't care about angle

% C0r= ((l^2+r^2)*(l^2+r0^2)) ^(1/4) /l;
C0z = (l^2+r0^2)^(1/4) /(l^2+r^2)^(1/4);
C1z= l/24*(  (3*r^2-2*l^2)/(l^2+r^2)^(3/2) +(2*l^2+9*r0^2)/(l^2+r0^2)^(3/2) ); 
% C1r= l/24*( (-2*l^2-9*r^2)/(l^2+r^2)^(3/2) +(2*l^2+9*r0^2)/(l^2+r0^2)^(3/2) ); 

pexi=r/r0*(l+sqrt(l^2+r0^2))/(l + sqrt(l^2+r^2) ) * exp(sqrt(l^2+r^2)/l) /exp(sqrt(l^2+r0^2)/l);
mexi=1/pexi;
% vr=0;
vz=0;
t=psih;
if(abs(r)<r0)
    tmp=1/( (mexi*exp(-1i*t))^nB -1 );
%     vr=            - 1/(2*pi*r)*C0r*imag(   tmp + C1r/nB*log(1+tmp) ); 
    vz= 1/(2*pi*l) + 1/(2*pi*l)*C0z*real(   tmp + C1z/nB*log(1+tmp) ); 
elseif(abs(r)>r0)
    tmp=1/( (pexi*exp(-1i*t))^nB -1 );
%     vr=            - 1/(2*pi*r)*C0r*imag(  tmp - C1r/nB*log(1+tmp) ); 
    vz=        0   + 1/(2*pi*l)*C0z*real( - tmp + C1z/nB*log(1+tmp) ); 
else
%     vr=0;
    vz=0;
end

% vt=l/r*( 1/(2*pi*l)-vz);
vz=sign*vz;
% ui=fact*nB*Gamma*[vr vt vz];
ui=fact*nB*Gamma*vz;
end


