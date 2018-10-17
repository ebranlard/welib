function [M K C]=getMatrices()
global Tower Shaft Rotor Nacelle Generator Algo Controller
pitch=Controller.pitch;
b=Rotor.Blade;
r=b.r;
ne=length(r);

%% Mass
%M_tot=Rotor.M+Nacelle.M;
%Mblade=Rotor.M/3;
Mblade=trapz(r,b.Mass);
M_tot=3*Mblade+Nacelle.M;   
% actually at one point the deflection of the blade should be taken into
% account instead of r


% Transformation matrices
[a12 a23 a34]=getTransfoMatrices(Controller.yaw, Nacelle.tilt, 0, Rotor.cone);
%pitch
a33b=[1 0 0;0 cosd(pitch) -sind(pitch);0 sind(pitch) cosd(pitch)];
a13b=a33b*a23*a12;
    

    
%%% modes vs modes
GM1=trapz(r,b.eigen1f(:,3).^2.*b.Mass + b.eigen1f(:,2).^2.*b.Mass  );
GM2=trapz(r,b.eigen1e(:,3).^2.*b.Mass + b.eigen1e(:,2).^2.*b.Mass  );
GM3=trapz(r,b.eigen2f(:,3).^2.*b.Mass + b.eigen2f(:,2).^2.*b.Mass  );

%%% modes vs x
xMass_in1=zeros(3,ne);
% horizontal rotor inertia
xMass_in1(3,:)=b.Mass;
%transformation of horizontal inertia in nacelle coordinates
xMass_in3b=(a13b*xMass_in1)';
IM1=trapz(r,xMass_in3b(:,3).*b.eigen1f(:,3) +xMass_in3b(:,2).*b.eigen1f(:,2) );
IM2=trapz(r,xMass_in3b(:,3).*b.eigen1e(:,3) +xMass_in3b(:,2).*b.eigen1e(:,2) );
IM3=trapz(r,xMass_in3b(:,3).*b.eigen2f(:,3) +xMass_in3b(:,2).*b.eigen2f(:,2) );

%%% modes vs theta
rotMass_in3=zeros(3,ne);
% rotational inertia without pitch
rotMass_in3(2,:)=b.Mass*cosd(Rotor.cone);
% rotational inertia with pitch
rotMass_in3b=(a33b*rotMass_in3)';
Rotor.I=3*trapz(r,b.Mass.*r.^2*cosd(Rotor.cone)^2);
RM1=trapz(r,r.*rotMass_in3b(:,3).*b.eigen1f(:,3)+r.*rotMass_in3b(:,2).*b.eigen1f(:,2) );
RM2=trapz(r,r.*rotMass_in3b(:,3).*b.eigen1e(:,3)+r.*rotMass_in3b(:,2).*b.eigen1e(:,2) );
RM3=trapz(r,r.*rotMass_in3b(:,3).*b.eigen2f(:,3)+r.*rotMass_in3b(:,2).*b.eigen2f(:,2) );

M(1,1)=M_tot;
M(2:3,2:3)=[Rotor.I+Generator.I Generator.I;Generator.I Generator.I];



%%%Code that works
%blade 1
% GM1=trapz([0;Rotor.r],[0;b.eigen1f(:,3).*b.Mass.*b.eigen1f(:,3)+b.eigen1f(:,2).*b.Mass.*b.eigen1f(:,2) ]  );
% GM2=trapz([0;Rotor.r],[0;b.eigen1e(:,3).*b.Mass.*b.eigen1e(:,3)+b.eigen1e(:,2).*b.Mass.*b.eigen1e(:,2) ]  );
% GM3=trapz([0;Rotor.r],[0;b.eigen2f(:,3).*b.Mass.*b.eigen2f(:,3)+b.eigen2f(:,2).*b.Mass.*b.eigen2f(:,2) ]  );
% 
% IM1=trapz([0;Rotor.r],[0;b.Mass.*b.eigen1f(:,3) ]  );
% IM2=trapz([0;Rotor.r],[0;b.Mass.*b.eigen1e(:,3) ]  );
% IM3=trapz([0;Rotor.r],[0;b.Mass.*b.eigen2f(:,3) ]  );
% 
% RM1=trapz([0;Rotor.r],[0;Rotor.r.*b.Mass.*b.eigen1f(:,2)]  );
% RM2=trapz([0;Rotor.r],[0;Rotor.r.*b.Mass.*b.eigen1e(:,2)]  );
% RM3=trapz([0;Rotor.r],[0;Rotor.r.*b.Mass.*b.eigen2f(:,2)]  );
% 
% M(1,1)=M_tot;
% M(2:3,2:3)=[Rotor.I+Generator.I Generator.I;Generator.I Generator.I];
%












% Blade / Blade
M(4:6,4:6)=[GM1 0  0;...
             0  GM2 0;...
             0   0  GM3];
M(7:9,7:9)=M(4:6,4:6);
M(10:12,10:12)=M(4:6,4:6);

% x / Blade
M(4:6,1)=[IM1;IM2;IM3];
M(7:9,1)=[IM1;IM2;IM3];
M(10:12,1)=[IM1;IM2;IM3];

% theta / Blade
M(4:6,2)=[RM1;RM2;RM3];
M(7:9,2)=[RM1;RM2;RM3];
M(10:12,2)=[RM1;RM2;RM3];

% nu / Blade
M(4:6,2)=[RM1;RM2;RM3];
M(7:9,2)=[RM1;RM2;RM3];
M(10:12,2)=[RM1;RM2;RM3];


% symmetry
M(1,:)=M(:,1);
M(2,:)=M(:,2);
M(3,:)=M(:,3);



%% stiffness matrix
K=zeros(12,12);
K(1,1)=Tower.k;
K(3,3)=Shaft.k;
Kblade=[ (2*pi*Rotor.Blade.eigenF(1))^2*GM1 0 0;...
         0 (2*pi*Rotor.Blade.eigenF(2))^2*GM2 0;...
         0 0 (2*pi*Rotor.Blade.eigenF(3))^2*GM3];
K(4:6,4:6)=Kblade;
K(7:9,7:9)=Kblade;
K(10:12,10:12)=Kblade;

%% Additional damping
C=zeros(12,12);
C(1,1)=10^3;
C(3,3)=Algo.damp;
delta=0.1;
Cblade=[ delta*2*Rotor.Blade.eigenF(1)*GM1 0 0;...
         0 delta*2*Rotor.Blade.eigenF(2)*GM2 0;...
         0 0 delta*2*Rotor.Blade.eigenF(3)*GM3];
C(4:6,4:6)=Cblade;
C(7:9,7:9)=Cblade;
C(10:12,10:12)=Cblade;


if(Algo.DontRotate)
    M(4:12,1:2)=0;
    M(1:2,4:12)=0;
end
end