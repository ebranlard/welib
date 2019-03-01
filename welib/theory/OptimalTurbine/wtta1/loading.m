%% Initialization / loading
clear all;
clc;
M=dlmread('airfoil.dat');
alpha_origin=M(:,1);
Cl_origin=M(:,2);
Cd_origin=M(:,3);
Profiles.Cl=M(:,2);
Profiles.Cd=M(:,3);
Profiles.alpha=M(:,1);
%% interpolation  on the range of interest
liminf=2;
limsup=32;
Cl=interp(Cl_origin(liminf:limsup),50);
Cd=interp(Cd_origin(liminf:limsup),50);
alpha=interp(alpha_origin(liminf:limsup),50);
Cl=[Cl_origin(1:(liminf-1)) ; Cl(1:(50*(limsup-liminf))) ; Cl_origin((limsup+1):33)];
Cd=[Cd_origin(1:(liminf-1)) ; Cd(1:(50*(limsup-liminf))) ; Cd_origin((limsup+1):33)];
alpha=[alpha_origin(1:(liminf-1)) ;  alpha(1:(50*(limsup-liminf))); alpha_origin((limsup+1):33)];

Profile.alpha=alpha;
Profile.Cl=Cl;
Profile.Cd=Cd;
ProfileNoDrag.alpha=alpha;
ProfileNoDrag.Cl=Cl;
ProfileNoDrag.Cd=Cd*0;

rho=1.225;
