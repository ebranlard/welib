function [ F ] = fGoldsteinOkulovF( H,B,r,N )
    R=1;
    w=1;
    Vr=(1/N:1/N:1)*R;
    KB=(Vr/R).^2./(H^2+(Vr/R).^2)  ;
    GammaGoldstein = fGoldtseinOkulovCirculation(H,w,R,B,Vr)*B/(H*w);
    % tip loss factor
    F=GammaGoldstein'./KB;
    F=interp1(Vr,F,r);
    
    
%     figure(12) 
%     hold all
%     plot(Vr/R,KB)
%     plot(Vr/R,GammaGoldstein)
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [ GammaGoldstein] = fGoldtseinOkulovCirculation(l,w,R,B,Vr)
%% Inititalization
l_bar=l/R;
VthetaB=(0:(B-1))*2*pi/B;  % [rad] Blade azimuthal angle
% Position of the helical vortices
N=length(Vr);
% Position of the Control Points
% The control points are placed between the vorticies
% They are used to compute velocities on the vortex sheet
VrCP=(3/(2*N):1/N:1)*R;   % [m] 

%% Calculation of matrices of influence at control points
% Circulation is unitary
%A_t=zeros(N,N); %tangential velocity
%A_z=zeros(N,N); %longitudinal velocity
A_x=zeros(N,N); %
for iB=1:B % loop on blades
    for j=1:length(Vr)  % loop on vorticies
        for i=1:length(VrCP) % loop on Control points
            a=Vr(j);  %[m] radius of helical vortex
            r=VrCP(i);%[m] radius of control point
%            A_t(i,j)=A_t(i,j) + fVt( r, a, l, -VthetaB(iB) )*1/(2*pi*l); %
%            A_z(i,j)=A_t(i,j) + fVz( r, a, l, -VthetaB(iB) )*1/(2*pi*l);
            A_x(i,j)=A_x(i,j) + fVx( r, a, l, -VthetaB(iB) )*1/(2*pi*l);
        end
    end
end
% Condition sum of gamma=0
%A_t(end,:)=1;
%A_z(end,:)=1;
A_x(end,:)=1;

% Calculation of boundary conditions values on the vortex sheet
% Boundary conditions are defined in equation 27, we apply then on the
% line made by the control points 
U_t0=zeros(N,1);
U_z0=zeros(N,1);
U_x0=zeros(N,1);
for i=1:(N-1)
     r=VrCP(i);      %[m]   
     r_bar=(r/R);   %[.]
     
     U_t0(i)= -w*r_bar*l_bar /(l_bar^2+r_bar^2 ); %l dimension less -> l_bar
     U_z0(i)=  w*r_bar^2     /(l_bar^2+r_bar^2 );
     U_x0(i)=  w*r/R /(l/R^2)*(1/(2*pi));
end
%% Solving for Gamma
%Gamma=A_z\U_z0; 
%Gamma2=A_t\U_t0; 
Gamma3=A_x\U_x0; 
GammaGoldstein=cumsum(Gamma3);

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
function [ Vt ] = fVt( r,a,l,t )
    %Vt=fVr(r,a,l,t)-r/l*fVz(r,a,l,t);
    Vt=l/r*(1-fVz(r,a,l,t));
end

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


