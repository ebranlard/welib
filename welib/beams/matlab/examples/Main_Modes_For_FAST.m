%% init 
clear all; clc; close all
addpath('../')
%%

PolyConstr=[1,1,1,1,1,0,0]; % Contraints for polyfit - NOTE: only 1s and then 0s for Modes.exe
nModes=3

% --- Tube L100-D8-t45- Diameter 8 -Thickness 45
D   = 8
t   = 0.045
L   = 100
rho = 7850                     ; % [kg/m3]
E   = 210e+09                  ; % [N/m^2] [Pa]  210 [GPa]
G   = 79.3e9                   ; % [N/m^2] [Pa]  Shear modulus. Steel
A   = pi*( (D/2)^2 - (D/2-t)^2);
Iy  = pi/64*(D^4-(D-2*t)^4)    ; % [m^4]
rho = 7850                     ; % [kg/m^3]
m   = rho*A                    ;
Ix  = pi/32*(D^4-(D-2*t)^4)    ; % Polar second moment of area [m^4]
Kt  = pi/64*(D^4-(D-2*t)^4)    ; % Torsion constant [m^4]
EI  = E*pi/64*(D^4-(D-2*t)^4)  ; % (Nm^2) I_annulus            = pi/4 (r2^4 - r1^4) = pi/64 (D2^4-D1^4)
Mtop=0;  

if Mtop>0
    [f,x,U,V,K] = fUniformBeamTheory('transverse-unloaded-topmass-clamped-free',EI,rho,A,L,'Mtop',Mtop,'norm','tip_norm');
else
    [f,x,U,V,K] = fUniformBeamTheory('transverse-unloaded-clamped-free',EI,rho,A,L,'norm','tip_norm');
end
x0=x/L; % IMPORTANT

for i=1:nModes
    fprintf('%.6f\t\n',f(i))
end



% --- Polyfit without constraints
p1=polyfit(x0,U(1,:),6);
p2=polyfit(x0,U(2,:),6);
% plot(x0, p(5)*x0.^2 + p(4)*x0.^3 + p(3)*x0.^4 + p(2)*x0.^5 + p(1)*x0.^6 ,'+' )
% plot(x0, p(5)*x0.^2 + p(4)*x0.^3 + p(3)*x0.^4 + p(2)*x0.^5 + p(1)*x0.^6 ,'+' )

% --- Polyfit with constraints
polyFunc = @(p) polyval(p.*PolyConstr,x0);
p1c = fminsearch( @(p)sum((U(1,:) - polyFunc(p)).^2), p1 );
p2c = fminsearch( @(p)sum((U(2,:) - polyFunc(p)).^2), p2 );

% --- Plot
figure,hold all
plot(x0,U(1,:))
plot(x0,U(2,:))
plot(x0,polyval(p1.*PolyConstr,x0),'+')
plot(x0,polyval(p2.*PolyConstr,x0),'+')
plot(x0,polyval(p1c.*PolyConstr,x0),'o')
plot(x0,polyval(p2c.*PolyConstr,x0),'o')
p1m=[-0.3263     1.2089 -1.4729 0.0042 1.5861  0.0 0.0]; % From Modes
p2m=[ 10.5859 -29.8978 20.0229 10.0191 -9.7301 0.0 0.0]; % From Modes
plot(x0,polyval(p1m.*PolyConstr,x0),'d')
plot(x0,polyval(p2m.*PolyConstr,x0),'d')


%% Modes input file
fprintf('False        Blade Switch:  True = Blade, False = Tower\n');
fprintf('     0.0     Steady state angular velocity of rotor (rpm)  [Ignored for towers]\n');
fprintf('     0.0     Pitch angle for blades (degrees)  [Ignored for towers]\n');
fprintf('%10.3f   Total beam length (m)\n',L);
fprintf('   0.0       Rigid beam length (m)\n');
fprintf('%10.1f   End mass (kg) \n',Mtop);
fprintf('%10d   Number of modes shapes or coefficients\n',sum(PolyConstr));
fprintf('%10d   Order of first coefficient\n',sum(PolyConstr==0));
fprintf('%10d   Number of input stations for distributed parameters\n',length(x0));
fprintf('   1.00      Factor to adjust mass\n');
fprintf('   1.00      Factor to adjust tower stiffness\n');
fprintf('   1.00      Factor to adjust in-plane stiffness  [Ignored for towers]\n');
for i= 1:length(x0)
    fprintf('%.5f %21.16e %21.13e\n',x0(i),m,EI)
end
fprintf('\n');
fprintf('Columns:\n')
fprintf('- Fractional distance along the flexible portion of the tower.  It must go from 0 to 1.\n')
fprintf('- mass/length in kg/m\n')
fprintf('- stiffness in N m^2.\n')
fprintf('---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------\n');
for i=2:6; fprintf('%13.6f     TwFAM1Sh(%d) - Mode 1, coefficient of x^%d term\n',p1c(7-i),i,i); end
for i=2:6; fprintf('%13.6f     TwFAM2Sh(%d) - Mode 2, coefficient of x^%d term\n',p2c(7-i),i,i); end
fprintf('-------------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------\n');
for i=2:6; fprintf('%13.6f     TwSSM1Sh(%d) - Mode 1, coefficient of x^%d term\n',p1c(7-i),i,i); end
for i=2:6; fprintf('%13.6f     TwSSM2Sh(%d) - Mode 2, coefficient of x^%d term\n',p2c(7-i),i,i); end

