function [ G ] = fTipLossGoldsteinOkulov( H,w,R,B,Vr,flagNoDim )
    %H is H okuilov, ie l_bar !!!!
    l=H*R;    
    G = fGoldtseinOkulovCirculation(l,w,R,B,Vr)*B/(H*w);
    if(flagNoDim)
        KB=fBetzCirculation(l,R,Vr);
        G = G./KB';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [ GammaGoldstein] = fGoldtseinOkulovCirculation(l,w,R,B,Vr)
drop=0;
if(Vr(1)==0)
   Vr=Vr(2:end);
   drop=1;
end
%% Inititalization
VthetaB=(0:(B-1))*2*pi/B;  % [rad] Blade azimuthal angle
% Position of the helical vortices
N=length(Vr);
% Position of the Control Points
% The control points are placed between the vorticies
% They are used to compute velocities on the vortex sheet
VrCP=(3/(2*N):1/N:1)*R;   % [m] 

%% Calculation of matrices of influence at control points for a helical band
% Circulation is unitary
A_x=zeros(N,N); %
for iB=1:B % loop on blades
    for j=1:length(Vr)  % loop on vorticies
        for i=1:length(VrCP) % loop on Control points
            a=Vr(j);  %[m] radius of helical vortex
            r=VrCP(i);%[m] radius of control point
            A_x(i,j)=A_x(i,j) + fVx( r, a, l, -VthetaB(iB) );%*1/(2*pi)*l^2;            
        end
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
     U_x0(i)= r/R*w;% w*r/R /((l/R)^2)*1/(2*pi);
end
%% Solving for Gamma
Gamma=A_x\U_x0; 
GammaGoldstein=cumsum(Gamma);
if(drop)
    GammaGoldstein=[0;GammaGoldstein];
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Vf ] = fVf(r,a,l,t  )
if(abs(r)<a)
    Vf=0 -l/r*(l^2+a^2)^(1/4) /(l^2+r^2)^(1/4) *( real( fSm(r/l,a/l,t) )+l/24*( (3*r^2-2*l^2)/(l^2+r^2)^(3/2) +(2*l^2+9*a^2)/(l^2+a^2)^(3/2) )*real(fSl(r/l,a/l,t)) ); 

else
    Vf=l/r -l/r*(l^2+a^2)^(1/4) /(l^2+r^2)^(1/4) *( -real( fSm(a/l,r/l,t) )+l/24*( (3*r^2-2*l^2)/(l^2+r^2)^(3/2) +(2*l^2+9*a^2)/(l^2+a^2)^(3/2) )*real(fSl(a/l,r/l,t)) ); 
end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ Vt ] = fVt( r,a,l,t )
%     %Vt=fVr(r,a,l,t)-r/l*fVz(r,a,l,t);
%     Vt=l/r*(1-fVz(r,a,l,t));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Vx ] = fVx( r,a,l,t )
    Vx=fVf(r,a,l,t)-r/l*fVz(r,a,l,t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vz] = fVz( r,a,l,t)
if(abs(r)<a)
    vz=1+ (l^2+a^2)^(1/4) /(l^2+r^2)^(1/4) *( real( fSm(r/l,a/l,t) )+l/24*( (3*r^2-2*l^2)/(l^2+r^2)^(3/2) +(2*l^2+9*a^2)/(l^2+a^2)^(3/2) )*real(fSl(r/l,a/l,t)) ); 
else
    vz=0+ (l^2+a^2)^(1/4) /(l^2+r^2)^(1/4) *( -real( fSm(a/l,r/l,t) )+l/24*( (3*r^2-2*l^2)/(l^2+r^2)^(3/2) +(2*l^2+9*a^2)/(l^2+a^2)^(3/2) )*real(fSl(a/l,r/l,t)) ); 
end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Sl ] = fSl( x,y,t )
Sl=-log(1-exp( (log(x/y * (sqrt(1+y^2)+1)/(sqrt(1+x^2)+1) )+sqrt(1+x^2)-sqrt(1+y^2))   +1i*t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Sm ] = fSm( x,y,t )
Sm=exp(1i*t)/( exp( -(  log(x/y * (sqrt(1+y^2) +1 )/(sqrt(1+x^2)+1)  ) +sqrt(1+x^2)-sqrt(1+y^2) ) )  - exp(1i*t)  );
end


