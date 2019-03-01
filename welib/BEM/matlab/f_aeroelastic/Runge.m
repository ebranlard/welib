
%% Initialization of time variables
dt=Algo.dt;
nt=floor(tmax/dt)+1;
ndof=Algo.DOF;  % degrees of freedom
NoUpdate=0;
Update=1;

t=zeros(1,nt);
% position, speed and acceleration
x=zeros(ndof,nt);
v=zeros(ndof,nt);
a=zeros(ndof,nt);

% Wind turbine loads
Torque=zeros(1,nt);
Thrust=zeros(1,nt);
Power=zeros(1,nt);
Pitch=zeros(1,nt);
HubV=zeros(1,nt);
MG=zeros(1,nt);
Flap=zeros(3,nt);
% initial position
x(:,1)=x0;
v(:,1)=v0;
a(:,1)=a0;

if(Algo.TwoDOF)
    x=x(1:2,:);
    v=v(1:2,:);
    a=a(1:2,:);
end
%% Time loop
INDEX=1:(nt-1);
RungeCore
%% errors
% mean((Thrust-Tower.k.*x(1,:)-M_tot*a(1,:))./Thrust)
% if(~Algo.TwoDOF)
%     mean((Torque-MG-a(2,:)*(Rotor.I+Generator.I)-a(3,:)*Generator.I   )./Torque)
% end
% mean((-MG-x(3,:)*Shaft.k-a(3,:)*Generator.I-a(2,:)*Generator.I    ))
