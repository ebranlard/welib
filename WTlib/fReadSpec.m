function [ WT ] = fReadSpec( WT, SpecFile )
fid=fopen(SpecFile);
A=textscan(fid,'%f %s %*[^\n]');
A=A{1};
fclose(fid);
WT.Rotor.nB=A(1);            %number of blades
WT.Rotor.BladeLength =A(2);   % BladeLength radius [m]
WT.Rotor.rhub = A(3);   % Hub  length along the coned line [m]
WT.Rotor.cone=A(4);  % cone angle [deg]
WT.Rotor.R =WT.Rotor.rhub+WT.Rotor.BladeLength;   %  Rotor unconed radius [m]
WT.Nacelle.tilt=A(5);
WT.Controller.yaw=A(6);
WT.Tower.H=A(7);   % height of the tower [m]
Algo.Ngrid = A(8);
 
WT.Spec.Omega_rated=A(9)*2*pi/60; % Nominal Rotational speed [rad/s]
WT.Spec.P_rated = A(10)/0.93; % kW

WT.Rotor.Omega=WT.Spec.Omega_rated;

end

