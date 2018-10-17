function [a12 a23 a34]=getTransfoMatrices(yaw,tilt,psi,cone)

a1=[ 1 0 0
    0 cosd(yaw) sind(yaw)
    0 -sind(yaw) cosd(yaw)];

a2=[cosd(tilt) 0 -sind(tilt)
    0 1 0
    sind(tilt) 0 cosd(tilt)];

a3=eye(3,3);

a12=a1*a2*a3;

a23=[cosd(psi) sind(psi) 0; 
    -sind(psi) cosd(psi) 0;
    0 0 1];

a34=[cosd(cone) 0 -sind(cone)
    0 1 0
    sind(cone) 0 cosd(cone)];
end