function P_trans=transformation2(theta_cone,theta_yaw,theta_tilt,theta_wing,P)
% Function transforms from the system 1 to the system 4.
% It is used to transfer grqavitational loads acting on the blades.
a1=[ 1          0                 0
     0 cos(theta_yaw) sin(theta_yaw)
     0 -sin(theta_yaw) cos(theta_yaw)];
a2=[cos(theta_tilt) 0 -sin(theta_tilt)
              0           1          0
    sin(theta_tilt) 0 cos(theta_tilt)];
a3=eye(3,3);
a12=a1*a2*a3;
a23=[cos(theta_wing) sin(theta_wing) 0    %transformation matrix to
transform from
    -sin(theta_wing) cos(theta_wing) 0    %tilt to shaft coordinate
system.
    0                    0            1];
a34=[cos(theta_cone) 0 -sin(theta_cone)
    0        1        0
    sin(theta_cone) 0 cos(theta_cone)];

P_trans=a34*a23*a12*P';
end
