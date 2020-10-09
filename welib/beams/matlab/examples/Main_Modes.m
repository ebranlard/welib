%% init 
clear all; clc; close all
addpath('../')
%%

L      = 100                      ;
EI     = 1.868E+12                ;
E      = 210e9                    ;
D      = 8                        ;
t      = 0.045                    ;
A      = pi*( (D/2)^2 - (D/2-t)^2);
EI     = E*pi/64*(D^4-(D-2*t)^4)  ; % (Nm^2) I_annulus= pi/4 (r2^4 - r1^4) = pi/64 (D2^4-D1^4)
rho    = 7850                     ;
M      = rho*A*L                  ; % Total mass of the beam
Pcr    = pi^2 *EI/(4*L^2)         ; % Euler critical buckling load for clamped free beam
omega0 = sqrt(EI/(rho*A*L^4))     ; % To make circular frequency dimensionless
f0     = omega0/(2*pi)            ; % to make frequency dimensionless
Mtop=1*M;
[f,x,U,V,K]     = fUniformBeamTheory('transverse-unloaded-clamped-free',EI,rho,A,L,'norm','tip_norm');%,varargin)
[fM,~,UM,VM,KM] = fUniformBeamTheory('transverse-unloaded-topmass-clamped-free',EI,rho,A,L,'Mtop',Mtop,'norm','tip_norm');%,varargin)

lambda=sqrt( 2*pi*fM*sqrt(rho*A/EI)) *L
figure,hold all
plot(x,U(1,:),x,UM(1,:))
plot(x,(1-cos(pi*x/(2*L))),'--');
%%
figure, hold all
plot(x,U(2,:),x,UM(2,:))
%%
for i=1:1
    fprintf('%.6f\t %.6f \n',f(i)/f0,fM(i)/f0)
end
