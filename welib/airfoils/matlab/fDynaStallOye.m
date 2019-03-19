function [Cl, fs]=fDynaStallOye(alpha_t,Polar,tau,fs_prev,dt)
% compute aerodynamical force from aerodynamic data
fs=0;
% dynamic stall
% interpolation from data
f_st  = interp1 ( Polar.alpha , Polar.f_st   , alpha_t )  ; 
Clinv = interp1 ( Polar.alpha , Polar.Cl_inv , alpha_t )  ; 
Clfs  = interp1 ( Polar.alpha , Polar.Cl_fs  , alpha_t )  ; 
% dynamic stall model

fs=f_st + (fs_prev-f_st )*exp(-dt/tau);


%Cl=fs*Clinv+(1-f_st)*Clfs; 
Cl=fs*Clinv+(1-fs)*Clfs; % N
% Cl  = interp1 ( Polar.alpha , Polar.Cl  , alpha_t )  ; 

