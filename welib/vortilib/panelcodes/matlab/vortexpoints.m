%% Initialization
clear all; close all; clc; %restoredefaultpath; %addpath(genpath())

%% Functions

function [U, V] = vortex_flow(DX, DY, Gamma, regParam, regMethod)
    % Induced velocity by one 2D vortex point on one Control Point (CP)
    % INPUTS:
    %   - DX, DY: vector from vortex to the control point (array-like)
    %   - Gamma: vortex strength
    %   - regParam: regularization parameter (default value: 0)
    %   - regMethod: regularization method (default value: [])
    % OUTPUTS:
    %   - U, V : Cartesian velocities (size of DX)
    if nargin < 4
        regParam = 0;
        regMethod = 0;
    end
    r2 = DX.^2 + DY.^2;
    tX = -DY;
    tY =  DX;
    bOK = r2 > 1e-8; % Avoid singularity at r=0
    U = zeros(size(DX));
    V = zeros(size(DY));
    if regMethod==0
        U(bOK) = Gamma/(2*pi) * tX(bOK)./r2(bOK); 
        V(bOK) = Gamma/(2*pi) * tY(bOK)./r2(bOK); 
    else
        U(bOK) = Gamma/(2*pi) * tX(bOK)./r2(bOK) .* (1 - exp(- r2(bOK) / regParam^2));
        V(bOK) = Gamma/(2*pi) * tY(bOK)./r2(bOK) .* (1 - exp(- r2(bOK) / regParam^2));
    end
end

function out = VP_panel_solve(x, y, Vinf_x, Vinf_y, hasLift)
    % Solve the flow about an airfoil using the vortex point method.
    % INPUTS:
    %  -x , y : airfoil points, assumed to go from TE to LE clockwise [m]
    %  -Vinf_xy: components of the freestream velocity in the x and y direction [m/s]
    % OUTPUTS:
    %  - out: storage for multiple variables, like Cp

    % Compute geometry (mid points, normals, tangents)
    P     = [x(:), y(:)]                   ;
    mids  = (P(1:end-1,:) + P(2:end,:)) / 2;
    dP    = P(2:end,:) - P(1:end-1,:)      ;
    ds    = sqrt(sum(dP.^2, 2))            ;
    t_hat = dP ./ ds                       ;
    n_hat = [-t_hat(:,2), t_hat(:,1)]      ;

    % Right hand side
    % We will implement the "Dirichlet" condition, no flow tangential "inside" of the body
    rhs = 0*x; % TODO Implement the right hand side: rhs = -Vinf . t
    rhs = -(Vinf_x*t_hat(:,1) + Vinf_y*t_hat(:,2));

    % Build matrix
    nCP = length(x)-1;  % Number of Control points
    nV = length(x)-1;   % Number of vortices
    CP = mids;          % Position of control points
    VP = mids;          % Position of vortices
    M = zeros(nCP,nV);
    for i = 1:nCP
        for j = 1:nV
            if i == j
                M(i,i) = 0.5/ds(j); % Exactly on panel, panel strength is half the velocity jump
            else
                ri = CP(i,:); % control point position
                rj = VP(j,:); % vortex point position
                % Induced velocity Vij from one unit vortex point at one control point.
                u = 0; % TODO Implement it
                v = 0; % TODO Implement it
                % Velocity projeted against the normal vector M[i,j] = Vij . t
                M(i,j) = 0; % TODO

                dr = ri - rj;
                [u, v] = vortex_flow(dr(1), dr(2), 1);
                % Velocity projected against the normal vector M(i,j) = Vij . t
                M(i,j) = (u*t_hat(i,1) + v*t_hat(i,2));
            end
            % Scale by length at the end (we formulate in terms of gammas)
            M(i,j) = M(i,j) * ds(j);
        end
    end
    if hasLift
        % Bonus: insert Kutta boundary condition here
    end

    % Solve (invert the system)
    gammas = M\rhs;
    if hasLift
        % Bonus: revert Kutta boundary condition here
    end
    Gammas = gammas.*ds;

    % Outputs
    out.VP = VP;            % Vortex points
    out.M = M;              % System matrix
    out.rhs = rhs;          % Right hand side
    out.Gammas = Gammas;    % Vortex points intensities

    % Velocity at wall (uses the fact that the vortex sheet creates a velocity jump)
    Vwall = zeros(size(mids));
    Vwall(:,1) = -(Gammas./ds).*t_hat(:,1);
    Vwall(:,2) = -(Gammas./ds).*t_hat(:,2);
    out.Vwall = Vwall;

    % Cp
    out.Cp = 0; % TODO Implement Cp using Vwall
    out.Cp = 1 - (Vwall(:,1).^2 + Vwall(:,2).^2) / (Vinf_x^2 + Vinf_y^2);
end



%% Main parameters
case_option = 'file';     % Choose between a file or the cylinder case
% case_option = 'cylinder'; % Choose between a file or the cylinder case
U0    = 1;                % Free stream [m/s]
alpha = 0*pi/180;         % Angle of attack [rad]

% Setup case
if strcmp(case_option, 'file')
    coords_filename = '../data/Diamond-coords.csv';
    data = xlsread(coords_filename);
    xa = data(:, 1); % Coordinates of airfoil
    ya = data(:, 2);
    Cp_theory = ones(length(xa)-1,1)*0.452
    hasLift = abs(alpha)>0;
elseif strcmp(case_option, 'cylinder')
    m         = 100;                % Number of panels
    R         = 2;                 % Cylinder radius [m]
    Gamma     = -4*pi*U0*R*0;      % Circulation - NOTE: Not available yet
    theta     = -linspace(0,2*pi, m+1) + asin(Gamma/(4*pi*U0*R));
    xa        = R*cos(theta);
    ya        = R*sin(theta);
    Cp_theory = 1-4*(sin((theta(1:end-1)+theta(2:end))/2) - Gamma/(4*pi*U0*R)).^2;
    hasLift = abs(Gamma)>0;
end

%% Derived parameters
Vinf_x = U0*cos(alpha);
Vinf_y = U0*sin(alpha);

%% Use the vortex panel method to find the solution
out = VP_panel_solve(xa, ya, Vinf_x, Vinf_y, hasLift);


if strcmp(case_option, 'file')
    fprintf('M matrix:\n')
    disp(out.M)
end


%% Plot results
figure; hold all;
plot(out.VP(:,1), Cp_theory, 'k-', 'DisplayName', 'Theory');
plot(out.VP(:,1), out.Cp, ':o', 'DisplayName', 'Numerical');
ylim([-4., 1.1]);
xlabel('x [m]');
ylabel('Cp [-]');
legend();
