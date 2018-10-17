%% Initialization of time variables
tplus=10;
tmax=tmax+tplus;
nt=length(t);
nplus=floor(tplus/dt)+1;

t=[t zeros(1,nplus)];
% position, speed and acceleration
x=[x zeros(ndof,nplus)];
v=[v zeros(ndof,nplus)];
a=[a zeros(ndof,nplus)];

% Wind turbine loads
Torque=[Torque zeros(1,nplus)];
Thrust=[Thrust zeros(1,nplus)];
Power=[Power zeros(1,nplus)];
MG=[MG zeros(1,nplus)];

%% Time loop
INDEX=nt:(nt+nplus-1);
RungeCore


%% errors
mean((Thrust-Tower.k.*x(1,:)-M_tot*a(1,:))./Thrust)
mean((Torque-MG-a(2,:)*(Rotor.I+Generator.I)-a(3,:)*Generator.I   )./Torque)
mean((-MG-x(3,:)*Shaft.k-a(3,:)*Generator.I-a(2,:)*Generator.I    ))