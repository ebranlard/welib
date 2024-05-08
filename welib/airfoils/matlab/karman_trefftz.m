% Set of static functions for the Karman Trefftz conformal map
classdef karman_trefftz
    methods (Static = true)
    function [z, dzdZ] = comf_map(Z, A, l)
        % Karman Trefftz conformal map from Z to z plane. Z: cylinder plane, z: airfoil plane
        z    = l*A*( (Z+A).^l + (Z-A).^l )./( (Z+A).^l - (Z-A).^l );
        dzdZ = 4*(l*A)^2 * (((Z-A).^(l-1)).*((Z+A).^(l-1)))./((((Z+A).^l)-((Z-A).^l)).^2);
    end 

    function [Z, dZdz] = comf_map_inv(z, A, l)
        % Inverse Karman Trefftz conformal map from z to Z plane. Z:cylinder plane, z:airfoil plane
        Z    = -A*(  ( ((z-l)./(z+l)).^(1/l) )+1  ) ./ ( (((z-l)./(z+l)).^(1/l))-1);
        dZdz = 1./((4*(l*A)^2)*(((Z-A).^(l-1)).*((Z+A).^(l-1)))./((((Z+A).^l)-((Z-A).^l)).^2));
    end 

    function [R, Beta, Gamma] =  cyl_params(XC, YC, A, U0, alpha)
        % Geometrical parameters of the Z-plane cylinder.
        R     = sqrt((A - XC)^2 + YC^2);
        Beta  = asin(- YC/R);
        Gamma = - 4 * pi * R * U0 * sin(alpha - Beta);
    end

    function [xa, ya] = shape(XC, YC, l, A, n)
        % Returns Karman-Trefftz profile coordinates
        % INPUT:
        %  - XC,YC: Center of Cylinder (in Z-plane)
        %  - l, A : Karman-Treffz parameters
        %  - n    : Number of points along mapped foil surface
        % OUTPUTS:
        %  - xa,ya: Coordinates of airfoil (in z-plane)

        % Cylinder parameters
        [R, Beta, Gamma] = karman_trefftz.cyl_params(XC, YC, A, 0, 0);
        % Cylinder coordinates in the Z-plane
        theta = linspace(Beta, -2*pi+Beta, n+1);
        theta = theta(1:n);
        ZC    = (XC + 1j * YC);
        Z_cyl = ZC + R * exp(1j * theta);
        % Transform back to the z-plane
        [za, dzdZ] = karman_trefftz.comf_map(Z_cyl, A, l);
        xa = real(za);
        ya = imag(za);
    end


    function [xa, ya, Cp_w, u_w, v_w] = wall(XC, YC, l, A, n, U0, alpha)
        % Returns the Karman-Trefftz profile coordinates and pressure coefficient
        % on the wall surface.
        %
        % INPUTS:
        %  - XC,YC: Center of Cylinder [m]
        %  - l, A : Karman-Treffz parameters
        %  - n    : Number of points along mapped surface
        %  - U0   : Freestream velocity [m/s]
        %  - alpha: angle of attack [rad]
        % OUTPUTS:
        %  - xa, ya:  coordinates of airfoil (in z-plane)
        %  - Cp_w : pressure coefficient at airfoil surface
        %  - u_w  : x-velocity at arifoil surface [m/s]
        %  - v_w  : y-velocity at arifoil surface [m/s]
        
        % Cylinder parameters
        [R, Beta, Gamma] =  karman_trefftz.cyl_params(XC, YC, A, U0, alpha);
        % Cylinder coordinates in the Z-plane
        theta = linspace(Beta, -2*pi+Beta, n+1);
        theta = theta(1:n);
        ZC    = (XC + 1j * YC);
        Z_cyl = ZC + R * exp(1j * theta);
        % Transform back to the z-plane
        [za, dzdZ] = karman_trefftz.comf_map(Z_cyl, A, l);
        xa = real(za);
        ya = imag(za);
        % Velocities on cylinder surface in Z-plane
        % TODO
        % Back in the z-plane
        % TODO
        u_w = xa*0
        v_w = xa*0
        Cp_w    = xa*0
    end

    function [u, v, Cp] = flow(x, y, XC, YC, l, A, U0, alpha)
        % Cylinder parameters
        R, Beta, Gamma =  karman_trefftz.cyl_params(XC, YC, A, U0, alpha)
        % Transform z to Z plane
        % TODO
        % Compute cylinder velocity in Z-plane
        % TODO
        % Back to z-plane
        % TODO
        u  = x*0 
        v  = x*0 
        Cp = x*0 
    end
end % methods
end % classdef
