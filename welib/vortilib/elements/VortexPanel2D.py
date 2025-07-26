"""

Constant Vortex Panel 2D (CVP)

References: 
 [1] Katz - Plotkin - Low speed aerodynamics p 270

"""

import numpy as np
import unittest
from scipy.integrate import quad, quad_vec
import math

# --------------------------------------------------------------------------------}
# --- Wrappers 
# --------------------------------------------------------------------------------{
def cvp_u(X, Y, SP1, SP2, gammas, debug=False):
    """ 
    Constant Vortex Panels, delimited by point 1 and 2 (can be discontinuous)

    INPUTS:
     - X: nd-array of control point x-coordinates
     - Y: nd-array of control point y-coordinates
     - SP1: (n,2)-array of first vortex panels points
     - SP2: (n,2)-array of second vortex panels points
     - gammas: (n)-array of vortex panels intensities (per unit length)
    """
    X = np.asarray(X)
    Y = np.asarray(Y)
    shp = X.shape
    X = X.flatten()
    Y = Y.flatten()
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    for j in range(len(SP1)):  # Loop on panels
        if debug:
            print('Panel', j, SP1[j], SP2[j], gammas[j])
        ui, vi = cvp_u1N(X, Y, SP1[j], SP2[j], gammas[j], debug=debug)
        U += ui
        V += vi
    U = U.reshape(shp)
    V = V.reshape(shp)
    return U, V

def ccvp_u(X, Y, SP, gammas, method=1, debug=False):
    """ 
    Contiguous Constant Vortex Panels 
    Contiguous => Panels are formed by consecutive points, can potentially form a closed loop

    INPUTS:
     - X: nd-array of control point x-coordinates
     - Y: nd-array of control point y-coordinates
     - SP: (n,2)-array of vortex panels points
     - gammas: (n-1)-array of vortex panels intensities (per unit length)
    """
    X = np.asarray(X)
    Y = np.asarray(Y)
    shp = X.shape
    X = X.flatten()
    Y = Y.flatten()
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    if False and method == 1:
       # Optimized for many control points
       for j in range(len(SP) - 1):  # Loop on panels
           ui, vi = cvp_u1N(X, Y, SP[j], SP[j + 1], gammas[j], debug=debug)
           U += ui
           V += vi
    else:
        # Not optimized, for any method
        for i in range(len(X)):
            for j in range(len(SP) - 1):  # Loop on panels
                ui, vi = cvp_u11((X[i], Y[i]), SP[j], SP[j + 1], gammas[j], method=method)
                U[i] += ui
                V[i] += vi
    U = U.reshape(shp)
    V = V.reshape(shp)
    return U, V

# --------------------------------------------------------------------------------}
# --- Vectorized Function for One Panel, Many Points
# --------------------------------------------------------------------------------{
def cvp_u1N(xCP, yCP, rS1, rS2, gamma=1, tol=1e-8, debug=False):
    """
    Velocity induced on N control points by one constant vortex panel (cvp)
    Formulae based on [1] p 273, but rotated in panel frame

    INPUTS:
     - xCP, yCP: position of control points, n-array
     - rS1: position of panel start (2 values)
     - rS2: position of panel end (2 values)
     - gamma: vortex strength per unit length
     - tol: tolerance for detecting if point is on panel
     - debug: if True, print debug information

    OUTPUTS:
     - U, V: velocities at control points, n-array
    """
    U = np.zeros_like(xCP)
    V = np.zeros_like(yCP)
    Un = np.zeros_like(xCP)
    Ut = np.zeros_like(yCP)
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    
    # --- Panel
    dx = xS2 - xS1
    dy = yS2 - yS1
    L = np.sqrt(dx**2 + dy**2)
    n_hat = np.array([-dy/L, dx/L])
    t_hat = np.array([ dx/L, dy/L])
    phi = np.atan2(dy, dx)
    phi = phi if phi>=0 else phi+2*np.pi

    # --- Vectors from panel points to CP
    x1C = xCP - xS1 
    y1C = yCP - yS1
    x2C = xCP - xS2
    y2C = yCP - yS2
    r1C_2 = x1C**2 + y1C**2
    r2C_2 = x2C**2 + y2C**2
    C = np.cos(phi)
    S = np.sin(phi)
    
    # ---Tilde coordinates
    # x_~ = x cos phi + y sin phi
    # y_~ =-x sin phi + y cos phi
    x_1C = x1C*C +  y1C*S  # xtilde_c - xtilde_1 = xtilde_c
    y_1C =-x1C*S +  y1C*C  # ytilde_c - ytilde_1 = ytilde_c
    x_2C = x2C*C +  y2C*S  # xtilde_c - xtilde_2
    y_2C =-x2C*S +  y2C*C  # = y_1C to machine precision

    # --- Angles measured in panel frame
    theta1 = np.arctan2(y_1C, x_1C)
    theta2 = np.arctan2(y_2C, x_2C)
    theta1[theta1 < 0] += 2 * np.pi
    theta2[theta2 < 0] += 2 * np.pi

    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot = dx * x1C + dy * y1C
    bOnLine = abs(cross) < tol
    bOnSeg = np.logical_and(-tol <= dot, dot <= L**2 + tol)
    b = np.logical_and(bOnLine, bOnSeg)

    # Principal value on the panel is gamma/2
    Un[b] = 0
    Ut[b] = gamma/2
    
    # Velocity off the panel (Katz & Plotkin), 10.21 - 10.22
    Un[~b] = gamma / (4 * np.pi) * (np.log(r2C_2[~b]/r1C_2[~b]))  # TODO sign
    Ut[~b] = gamma / (2 * np.pi) * (theta2[~b] - theta1[~b])

    U = Un*n_hat[0] + Ut*t_hat[0]
    V = Un*n_hat[1] + Ut*t_hat[1]
    
    if debug:
        print('xCP ', xCP, yCP)
        print('xS1 ', xS1, yS1)
        print('x1C', x1C, y1C)
        print('L, dx, dy', L, dx, dy)
        print('cross, dot', cross, dot)
        print('t_hat', t_hat, gamma, b)
    return U, V

# --------------------------------------------------------------------------------}
# --- 11 : Velocity (u) induced by one (1) constant vortex panel (cvp) on one control point (1)
# --------------------------------------------------------------------------------{
def cvp_u11(rCP, rS1, rS2, gamma=1, method=1, tol=1e-8, principal=False, WARN=[0]):
    """
    Velocity (u) induced by one (1) constant vortex panel (cvp) on one control point (1)

    Global frame: X, Y, panel frame x, y
      x = X cos phi + Y sin phi
      y =-X sin phi + Y cos phi

    INPUTS:
    - rCP: position of control point (two values)
    - rS1: position of panel start (two values)
    - rS2: position of panel end (two values)
    - gamma: vortex strength per unit length
    - method: 0 (Anderson), 1 (Theoretical), 2 (Numerical Quadrature), 10 (Point vortex sum)
    - principal: if True, return principal value no matter what (if we know we are on a panel)
                 Otherwise, the code will attempt to detect if we are on the panel.

    OUTPUTS:
    - u, v: velocity components at control point
    """
    if method == 0:
        # TODO principal value not implemented
        u, v = cvp_u11_anderson(rCP, rS1, rS2, gamma=gamma, tol=tol)
    elif method == 1:
        u, v = cvp_u11_kp(rCP, rS1, rS2, gamma=gamma, tol=tol, principal=principal)
    elif method == 2:
        u, v = cvp_u11_quad(rCP, rS1, rS2, gamma=gamma, tol=tol, principal=principal)
    elif method == 10:
        from welib.vortilib.elements.VortexPoint import vp_u
        # TODO principal value not implemented
        # We use many point vortices along the panel
        nS = 30
        xS1, yS1 = rS1
        xS2, yS2 = rS2
        xCP, yCP = rCP
        dx = xS2 - xS1
        dy = yS2 - yS1
        L = np.sqrt(dx ** 2 + dy ** 2)
        dl = L / nS
        Gamma = gamma * dl
        P = np.linspace(rS1, rS2, nS + 1)
        Ps = (P[:-1, :] + P[1:, :]) / 2
        u, v = 0, 0
        for i in range(nS):
            ui, vi = vp_u(xCP, yCP, (Ps[i][0], Ps[i][1]), Gamma=Gamma)
            u += ui
            v += vi
    else:
        raise ValueError("Method must be 0 (Anderson), 1 (Theoretical), 2 (Quadrature) or 10 (Rieman sum)")

    return u, v

def cvp_u11_quad(rCP, rS1, rS2, gamma=1, tol=1e-8, principal=False):
    """
    Velocity induced by a constant vortex panel using numerical quadrature.

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end (two values)
     - gamma: vortex strength per unit length
     - tol: tolerance for detecting if point is on panel
     - principal: if True, return principal value

    OUTPUTS:
     - u, v: velocity components at control point
    """
    xCP, yCP = rCP
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    dx = xS2 - xS1
    dy = yS2 - yS1
    L = np.sqrt(dx**2 + dy**2)
    n_hat = np.array([-dy/L, dx/L])
    t_hat = np.array([ dx/L, dy/L])
    
    # Check if point is on the panel
    dX_r1 = xCP - xS1
    dY_r1 = yCP - yS1
    cross = dx * dY_r1 - dy * dX_r1
    dot = dx * dX_r1 + dy * dY_r1
    
    phi = np.atan2(dy, dx)
    phi = phi if phi >= 0 else phi + 2 * np.pi
    
    # --- Vectors from panel points to CP
    x1C = xCP - xS1
    y1C = yCP - yS1
    x2C = xCP - xS2
    y2C = yCP - yS2
    r1C_2 = x1C**2 + y1C**2
    r2C_2 = x2C**2 + y2C**2
    C = np.cos(phi)
    S = np.sin(phi)
    
    # ---Tilde coordinates
    # x_~ = x cos phi + y sin phi
    # y_~ =-x sin phi + y cos phi
    x_1C = x1C*C +  y1C*S  
    y_1C =-x1C*S +  y1C*C
    x_2C = x2C*C +  y2C*S
    y_2C =-x2C*S +  y2C*C

    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot   = dx * x1C + dy * y1C
    if principal or (abs(cross) < tol and 0 - tol <= dot <= L**2 + tol):
        ut = gamma/2
        un = 0
        u, v = un * n_hat + ut * t_hat 
        return u, v
    
    def integrand_un_ut(s):
        ss = s/L
        x_minus_x0 = x_1C*(1-ss) + x_2C * ss
        denom =  (x_minus_x0**2 + y_1C**2 )
        ut = - gamma / (2 * np.pi) * y_1C / denom
        un =   gamma / (2 * np.pi) * x_minus_x0 / denom # TODO sign
        return np.array( [un, ut] )
    result, _ = quad_vec(integrand_un_ut, 0, L)
    un, ut = result[0], result[1]
    u, v = un * n_hat + ut * t_hat 
    return u, v

def cvp_u11_kp(rCP, rS1, rS2, gamma=1, tol=1e-8, principal=False):
    """
    Velocity induced on 1 control point by one vortex panel

    Formulae based on [1] p 272, but rotated in panel frame
    
    Global frame: x, y, panel frame x~, y~ (tilde)
      x~ = x cos phi + y sin phi
      y~ =-x sin phi + y cos phi

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end (two values)
     - gamma: vortex strength per unit length
     - tol: tolerance for detecting if point is on panel
     - principal: if True, return principal value (gamma/2 tangential)

    OUTPUTS:
     - u, v: velocity components at control point
    """
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    xCP, yCP = rCP
    
    # --- Panel
    dx = xS2 - xS1
    dy = yS2 - yS1
    L = np.sqrt(dx**2 + dy**2)
    n_hat = np.array([-dy/L, dx/L])
    t_hat = np.array([ dx/L, dy/L])
    phi = np.atan2(dy, dx)
    phi = phi if phi >= 0 else phi + 2 * np.pi
    
    # --- Vectors from panel points to CP
    x1C = xCP - xS1
    y1C = yCP - yS1
    x2C = xCP - xS2
    y2C = yCP - yS2
    r1C_2 = x1C**2 + y1C**2
    r2C_2 = x2C**2 + y2C**2
    C = np.cos(phi)
    S = np.sin(phi)
    
    # ---Tilde coordinates
    # x_~ = x cos phi + y sin phi
    # y_~ =-x sin phi + y cos phi
    x_1C = x1C*C +  y1C*S  # xtilde_c - xtilde_1 = xtilde_c
    y_1C =-x1C*S +  y1C*C  # ytilde_c - ytilde_1 = ytilde_c
    x_2C = x2C*C +  y2C*S  # xtilde_c - xtilde_2
    y_2C =-x2C*S +  y2C*C  # = y_1C to machine precision

    # --- Angles measured in panel frame    
    theta1 = np.arctan2(y_1C, x_1C)
    theta2 = np.arctan2(y_2C, x_2C)
    if theta2<0: 
        theta2 += 2*np.pi # Panel angles are assumed positive
    if theta1<0: 
        theta1 += 2*np.pi # Panel angles are assumed positive
    
    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot   = dx * x1C + dy * y1C
    if principal or (abs(cross) < tol and 0 - tol <= dot <= L**2 + tol):
        ut = gamma/2
        un = 0
        u,v = un*n_hat + ut*t_hat
        return u, v
    # Velocity off the panel (Katz & Plotkin), 10.39 - 10.40
    un =  - gamma / (4 * np.pi) *(np.log(r2C_2/r1C_2)) # NOTE: using opposive convention for sign
    ut =  - gamma / (2 * np.pi) * (theta2 - theta1)
    u,v = un*n_hat + ut*t_hat
    return u, v

def cvp_u11_anderson(rCP, rS1, rS2, gamma=1, tol=1e-8, principal=False):
    """
    Velocity induced by a single constant-strength vortex panel using Anderson's panel influence formula.

    INPUTS:
    - rCP: position of control point (2 values)
    - rS1: position of panel start (2 values)
    - rS2: position of panel end (2 values)
    - gamma: vortex strength per unit length

    OUTPUTS:
    - u, v: velocity components at control point
    """
    xCP, yCP = rCP
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    dx       = xS2 - xS1
    dy       = yS2 - yS1
    L        = np.sqrt(dx**2 + dy**2)
    phi      = np.atan2(dy, dx)
    cos_phi  = np.cos(phi)
    sin_phi  = np.sin(phi)

    # Panel start point
    x1C = xCP - xS1
    y1C = yCP - yS1

    # Check if the point is on the panel (principal value)
    cross = dx * y1C - dy * x1C
    dot   = dx * x1C + dy * y1C
    if principal or (abs(cross) < tol and 0 - tol <= dot <= L**2 + tol):
        # Tangential and normal unit vectors
        n_hat = np.array([-dy/L, dx/L])
        t_hat = np.array([ dx/L, dy/L])
        ut    = gamma/2
        un    = 0
        u, v  = un*n_hat + ut*t_hat
        return u, v

    # Anderson's notation
    A  = -x1C * cos_phi - y1C * sin_phi
    B  = x1C**2 + y1C**2
    Cx = sin_phi
    Dx = -y1C
    Cy = -cos_phi
    Dy = x1C
    E  = np.sqrt(max(B - A**2, 0.0))

    if E < tol or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
        Nx = 0.0
        Ny = 0.0
    else:
        denom = L**2 + 2*A*L + B
        if denom <= 0 or B <= 0:
            Nx = 0.0
            Ny = 0.0
        else:
            atan_diff = np.arctan2(L + A, E) - np.arctan2(A, E)
            Nx = 0.5 * Cx * np.log(denom / B) + ((Dx - A*Cx) / E) * atan_diff
            Ny = 0.5 * Cy * np.log(denom / B) + ((Dy - A*Cy) / E) * atan_diff
            if np.iscomplex(Nx) or np.isnan(Nx) or np.isinf(Nx):
                Nx = 0.0
            if np.iscomplex(Ny) or np.isnan(Ny) or np.isinf(Ny):
                Ny = 0.0
    u = gamma / (2 * np.pi) * Nx
    v = gamma / (2 * np.pi) * Ny
    return u, v

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_CVP_flow(self, plot=False):
        """
        Test velocity field induced by a vortex panel using different methods and optionally plot the results.
        """
        from welib.CFD.flows2D import flowfield2D, flowfield2D_plot
        gamma = 1
        PP1 = np.array([-0.5, -0.3])
        PP2 = np.array([0.5, 0.4])
        SP = np.vstack((PP1, PP2))

        xmax = 1.5
        xs = np.linspace(-xmax, xmax, 15)
        dy = 0.5
        ys = xs * 0 + dy
        xs = np.concatenate((xs, xs))
        ys = np.concatenate((ys, ys - 2 * dy))

        vel0 = lambda X, Y: ccvp_u(X, Y, SP, [gamma], method=0)
        vel1 = lambda X, Y: ccvp_u(X, Y, SP, [gamma], method=1)
        vel2 = lambda X, Y: ccvp_u(X, Y, SP, [gamma], method=2)
        vel4 = lambda X, Y: ccvp_u(X, Y, SP, [gamma], method=10)

        X, Y, U0, V0 = flowfield2D(vel0, xmax=1.5, ymin=-1.3, nx=15)
        X, Y, U1, V1 = flowfield2D(vel1, xmax=1.5, ymin=-1.3, nx=15)
        X, Y, U2, V2 = flowfield2D(vel2, xmax=1.5, ymin=-1.3, nx=15)
        X, Y, U4, V4 = flowfield2D(vel4, xmax=1.5, ymin=-1.3, nx=15)

        np.testing.assert_almost_equal(U1, U2, decimal=6)
        np.testing.assert_almost_equal(V1, V2, decimal=6)
        #np.testing.assert_almost_equal(U1, U4, decimal=1)
        #np.testing.assert_almost_equal(V1, V4, decimal=1)

        if plot:
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 4, sharey=False, figsize=(10, 3.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.50)
            ax =  flowfield2D_plot(X, Y, U1, V1, bounded=False, xs=xs, ys = ys, ax=axes[0], maxVal=0.5, minVal=0)
            ax.set_title('Theoretical')

            ax = flowfield2D_plot(X, Y, U0, V0, bounded=False, xs=xs, ys=ys, ax=axes[1], maxVal=0.5, minVal=0)
            ax.set_title('Anderson')

            ax = flowfield2D_plot(X, Y, U2, V2, bounded=False, xs=xs, ys=ys, ax=axes[2], maxVal=0.5, minVal=0)
            ax.set_title('Quadrature')

            ax = flowfield2D_plot(X, Y, U4, V4, bounded=False, xs=xs, ys=ys, ax=axes[3], maxVal=0.5, minVal=0)
            ax.set_title('Point sum')
            plt.show()

    def test_CVP_PrincipalValue(self):
        # --- One Panel on the x-axis - Check that principal value is obtained
        gamma = 10
        PP1 = [-0.5, 0.0]
        PP2 = [0.5, 0.0]
        rCP = [0., 0]
        L = 1

        u0 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=0)
        u1 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=1)
        u2 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=2)
        u4 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=10)
        np.testing.assert_almost_equal(u0[0], gamma/2, decimal=6)
        np.testing.assert_almost_equal(u1[0], gamma/2, decimal=6)
        np.testing.assert_almost_equal(u2[0], gamma/2, decimal=6)
        #np.testing.assert_almost_equal(u4[0], gamma/2, decimal=6) # No hope
        np.testing.assert_almost_equal(u0[1], 0, decimal=6)
        np.testing.assert_almost_equal(u1[1], 0, decimal=6)
        np.testing.assert_almost_equal(u2[1], 0, decimal=6)
        #np.testing.assert_almost_equal(u4[1], 0, decimal=6)

        # --- One panel on x-axis - principal value with many points
        SP = np.vstack((PP1, PP2))
        x = np.linspace(-0.25, 0.25, 10)
        y = x * 0 + 0
        u0 = ccvp_u(x, y, SP, gammas=[gamma], method=0)
        u1 = ccvp_u(x, y, SP, gammas=[gamma], method=1)
        u2 = ccvp_u(x, y, SP, gammas=[gamma], method=2)
        #u4 = ccvp_u(x, y, SP, gammas=[gamma], method=10)
        u_ref = ([gamma/2]*len(x), [0] * len(x))
        np.testing.assert_almost_equal(u0, u_ref, decimal=6)
        np.testing.assert_almost_equal(u1, u_ref, decimal=6)
        np.testing.assert_almost_equal(u2, u_ref, decimal=6)
        #np.testing.assert_almost_equal(u4, u_ref, decimal=6)
        
        # --- One tilted panel - principal value with one point
        PP1 = [-0.5, 0.2]
        PP2 = [0.5, 0.4]
        rCP = [(PP1[0] + PP2[0])/2, (PP1[1] + PP2[1])/2]  # Midpoint
        u0 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=0, principal=True)
        u1 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=1, principal=True)
        u2 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=2, principal=True)
        #u4 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=10, principal=True)
        np.testing.assert_almost_equal(u1, u2, decimal=6)
        np.testing.assert_almost_equal(u1, u0, decimal=6)
        #np.testing.assert_almost_equal(u1, u4, decimal=6)
    def test_CVP_crossing_points(self, plot=False):
        """
        Test velocity induced by a tilted panel at points crossing it at an angle through its midpoint.
        Compares methods 1 (Theoretical) and 2 (Quadrature).
        """
        gamma = 3
        PP1 = np.array([-0.5, 0.2])
        PP2 = np.array([0.5, 0.4])
        #PP1 = np.array([-0.5, 0.0])
        #PP2 = np.array([0.5, 0.0])
        SP = np.vstack((PP1, PP2))
        dx = PP2[0] - PP1[0]
        dy = PP2[1] - PP1[1]
        L = np.sqrt(dx**2 + dy**2)
        midpoint = np.array([(PP1[0] + PP2[0])/2, (PP1[1] + PP2[1])/2])

        PV_mid = gamma/2

        # Define points crossing the panel at 45 degrees through the midpoint
        n_points = 41
        s = np.linspace(-1, 1, n_points)
        phi = np.atan2(dy, dx)
        for angle_off in [0, np.pi/4, np.pi/2]:
            angle = phi + angle_off
            x = midpoint[0] + s * np.cos(angle)
            y = midpoint[1] + s * np.sin(angle)
            rCPs = np.vstack((x, y)).T

            # Compute velocities for each method
            u0, v0 = np.zeros(n_points), np.zeros(n_points)
            u1, v1 = np.zeros(n_points), np.zeros(n_points)
            u2, v2 = np.zeros(n_points), np.zeros(n_points)
            u4, v4 = np.zeros(n_points), np.zeros(n_points)
            for i, rCP in enumerate(rCPs):
                u0[i], v0[i] = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=0)
                u1[i], v1[i] = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=1)
                u2[i], v2[i] = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=2)
                u4[i], v4[i] = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=10)

            # Verify principal value at midpoint (s=0, index=n_points//2)
            expected = (np.cos(phi)*PV_mid, np.sin(phi)*PV_mid)
            np.testing.assert_almost_equal([u0[n_points//2], v0[n_points//2]], expected, decimal=6)
            np.testing.assert_almost_equal([u1[n_points//2], v1[n_points//2]], expected, decimal=6)
            np.testing.assert_almost_equal([u2[n_points//2], v2[n_points//2]], expected, decimal=6)
            #np.testing.assert_almost_equal([u4[n_points//2], v4[n_points//2]], expected, decimal=6) # No hope

            # Compare methods
            np.testing.assert_almost_equal(u1, u2, decimal=6)
            #np.testing.assert_almost_equal(u1, u0, decimal=6)
            #np.testing.assert_almost_equal(u1, u4, decimal=6) # can't get principal value
            #np.testing.assert_almost_equal(v1, v0, decimal=6)
            np.testing.assert_almost_equal(v1, v2, decimal=6)
            #np.testing.assert_almost_equal(v1, v4, decimal=6) # can't get principal value

            if plot:
                import matplotlib.pyplot as plt
                # Plot 1: Panel and control points
                fig, ax = plt.subplots(figsize=(6, 6))
                ax.plot([PP1[0], PP2[0]], [PP1[1], PP2[1]], 'b-', label='Panel')
                ax.scatter(x, y, c='r', s=50, label='Control Points')
                ax.scatter(midpoint[0], midpoint[1], c='g', s=100, marker='x', label='Midpoint')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_title('Panel and Crossing Points')
                ax.legend()
                ax.axis('equal')

                # Plot 2: u and v velocities
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
                ax1.plot(s, u0, 'k.', label='Anderson')
                ax1.plot(s, u1, 'r--', label='Theoretical')
                ax1.plot(s, u2, 'g:', label='Quadrature')
                ax1.plot(s, u4, 'b-.', label='Riemann sum')
                ax1.set_ylabel('u velocity')
                ax1.legend()
                ax1.grid(True)
                ax2.plot(s, v0, 'k.', label='Anderson')
                ax2.plot(s, v1, 'r--', label='Theoretical')
                ax2.plot(s, v2, 'g:', label='Quadrature')
                ax2.plot(s, v4, 'b-.', label='Rieman sum')
                ax2.set_xlabel('s (along crossing line)')
                ax2.set_ylabel('v velocity')
                ax2.legend()
                ax2.grid(True)
                plt.tight_layout()
        if plot:
            plt.show()

    def test_CVP_flowrate(self):
        # Test flow rate and circulation for a vortex panel
        from welib.CFD.flows2D import flowrate2D
        from welib.CFD.flows2D import circulation2D

        gamma = 3
        PP1 = np.array([-0.6, 0.0])
        PP2 = np.array([ 0.6, 0.0])
        L = 1.2 
        SP = np.vstack((PP1, PP2))
        theta = np.linspace(0, 2 * np.pi, 500, endpoint=True)
        # Contour over all of the panel, we should get 0
        for R in [1.5]:
            x = R * np.cos(theta)
            y = R * np.sin(theta)
            u0, v0 = ccvp_u(x, y, SP, gammas=[gamma], method=0)
            u1, v1 = ccvp_u(x, y, SP, gammas=[gamma], method=1)
            u2, v2 = ccvp_u(x, y, SP, gammas=[gamma], method=2)
            u4, v4 = ccvp_u(x, y, SP, gammas=[gamma], method=10)

            Gamma0 = circulation2D(x, y, u0, v0, verbose=False)
            Gamma1 = circulation2D(x, y, u1, v1, verbose=False)
            Gamma2 = circulation2D(x, y, u2, v2, verbose=False)
            Gamma4 = circulation2D(x, y, u4, v4, verbose=False)
            Q0 = flowrate2D(x, y, u0, v0, verbose=False, ns=-1)
            Q1 = flowrate2D(x, y, u1, v1, verbose=False, ns=-1)
            Q2 = flowrate2D(x, y, u2, v2, verbose=False, ns=-1)
            Q4 = flowrate2D(x, y, u4, v4, verbose=False, ns=-1)
            np.testing.assert_almost_equal(Gamma0, gamma*L, decimal=2) # <<<
            np.testing.assert_almost_equal(Gamma1, gamma*L, decimal=4)
            np.testing.assert_almost_equal(Gamma2, gamma*L, decimal=4)
            np.testing.assert_almost_equal(Gamma4, gamma*L, decimal=4)
            np.testing.assert_almost_equal(Q0, 0, decimal=8)
            np.testing.assert_almost_equal(Q1, 0, decimal=8)
            np.testing.assert_almost_equal(Q2, 0, decimal=8)
            np.testing.assert_almost_equal(Q4, 0, decimal=8)

    def test_CVP_debug_point(self):
        """
        Debug test to compare velocities at a specific off-panel point.
        """
        gamma = 1
        PP1 = np.array([-0.5, 0.0])
        PP2 = np.array([0.5, 0.0])
        rCP = [0.0, 0.5]  # Point above the panel
        u0, v0 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=0)
        u1, v1 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=1)
        u2, v2 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=2)
        u4, v4 = cvp_u11(rCP, PP1, PP2, gamma=gamma, method=10)
        #print(f"Method 1  (Theoretical): u={u1}, v={v1}")
        #print(f"Method 2  (Quadrature) : u={u2}, v={v2}")
        #print(f"Method 0  (Anderson  ) : u={u0}, v={v0}")
        #print(f"Method 10 (Rieman sum) : u={u4}, v={v4}")
        np.testing.assert_almost_equal([u1, v1], [u2, v2], decimal=6)
        np.testing.assert_almost_equal([u1, v1], [u0, v0], decimal=6)
        np.testing.assert_almost_equal([u1, v1], [u4, v4], decimal=3)


if __name__ == "__main__":
    #Test().test_CVP_debug_point()
    #Test().test_CVP_flow(plot=True)
    #Test().test_CVP_PrincipalValue()
    #Test().test_CVP_crossing_points(plot=True)
    #Test().test_CVP_flowrate()
    unittest.main()
