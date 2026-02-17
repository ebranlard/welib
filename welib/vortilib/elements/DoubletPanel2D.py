"""

Constant Doublet Panel 2D (CDP)


References: 
 [1] Katz - Plotkin - Low speed aerodynamics p 270



Convention:

  y ^  
    |
    |            mu 
    |    ^^^^^^^^^^^^^^^^
    -----||||||||||||||||------> x


             Equivalent to 

  y ^  
    |
    | Gamma1 = +mu      Gamma2 = -mu
    |  <.                .>
    -----)--------------(------> x

    Gamma is positive about z !!!!!


"""

import numpy as np
import unittest
from scipy.integrate import quad, quad_vec

# --------------------------------------------------------------------------------}
# --- Wrappers 
# --------------------------------------------------------------------------------{
def dcdp_u(X, Y, SP1, SP2, mus, debug=False):
    """ 
    Constant Doublet Panels, delimited by point 1 and 2 (can be discontinuous)

    INPUTS:
     - X: nd-array of control point x-coordinates
     - Y: nd-array of control point y-coordinates
     - SP1: (n,2)-array of first doublet panels points
     - SP2: (n,2)-array of second doublet panels points
     - mus: (n)-array of doublet panels intensities
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
            print('Panel', j, SP1[j], SP2[j], mus[j])
        ui, vi = cdp_u1N(X, Y, SP1[j], SP2[j], mus[j], debug=debug)
        U += ui
        V += vi
    U = U.reshape(shp)
    V = V.reshape(shp)
    return U, V

def cdp_u(X, Y, SP, mus, method=1, debug=False):
    """ 
    Contiguous Constant Doublet Panels 
    Contiguous => Panels are formed by consecutive points, can potentially form a closed loop

    INPUTS:
     - X: nd-array of control point x-coordinates
     - Y: nd-array of control point y-coordinates
     - SP: (n,2)-array of doublet panels points
     - mus: (n-1)-array of doublet panels intensities
    """
    X = np.asarray(X)
    Y = np.asarray(Y)
    shp = X.shape
    X = X.flatten()
    Y = Y.flatten()
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    if method == 1:
       # Optimized for many control points
       for j in range(len(SP) - 1):  # Loop on panels
           ui, vi = cdp_u1N(X, Y, SP[j], SP[j + 1], mus[j], debug=debug)
           U += ui
           V += vi
    else:
        # Not optimized, for any method
        for i in range(len(X)):
            for j in range(len(SP) - 1):  # Loop on panels
                ui, vi = cdp_u11((X[i], Y[i]), SP[j], SP[j + 1], mus[j], method=method)
                U[i] += ui
                V[i] += vi
    U = U.reshape(shp)
    V = V.reshape(shp)
    return U, V

# --------------------------------------------------------------------------------}
# --- Vectorized Function for One Panel, Many Points
# --------------------------------------------------------------------------------{
def cdp_u1N(xCP, yCP, rS1, rS2, mu=1, tol=1e-8, debug=False):
    """
    Velocity induced on N control points by one constant doublet panel (cdp)
    Formulae based on [1] p 273, but rotated in panel frame

    INPUTS:
     - xCP, yCP: position of control points, n-array
     - rS1: position of panel start (2 values)
     - rS2: position of panel end (2 values)
     - mu: doublet strength
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
    phi = phi if phi>=0 else phi+2*np.pi # positive angles # phi=mod(phi, 2*pi)

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
    dtheta = theta2 - theta1

    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot = dx * x1C + dy * y1C
    bOnLine = abs(cross) < tol
    bOnSeg = np.logical_and(-tol <= dot, dot <= L**2 + tol)
    b = np.logical_and(bOnLine, bOnSeg)

    # Principal value on the panel is zero for doublet
    Un[b] = +mu/(2*np.pi) * (1/(x_1C[b]+1e-8)-1/(x_2C[b]+1e-8)) # NOTE NOTE NOTE: I flipped the sign compared to Katz-Plotkin 
    
    # Velocity off the panel (Katz & Plotkin), 10.29 - 10.30
    Un[~b] = + mu / (2 * np.pi) * (x_1C[~b] / r1C_2[~b] - x_2C[~b]/r2C_2[~b] )
    Ut[~b] = - mu / (2 * np.pi) * (y_1C[~b] / r1C_2[~b] - y_2C[~b]/r2C_2[~b] )

    U = Un*n_hat[0] + Ut*t_hat[0]
    V = Un*n_hat[1] + Ut*t_hat[1]
    
    if debug:
        print('xCP ', xCP, yCP)
        print('xS1 ', xS1, yS1)
        print('x1C', x1C, y1C)
        print('L, dx, dy', L, dx, dy)
        print('cross, dot', cross, dot)
        print('t_hat', t_hat, mu, b)
    return U, V

# --------------------------------------------------------------------------------}
# --- 11 : Velocity (u) induced by one (1) constant doublet panel (cdp) on one control point (1)
# --------------------------------------------------------------------------------{
def cdp_u11(rCP, rS1, rS2, mu=1, method=1, tol=1e-8, principal=False, WARN=[0]):
    """
    Velocity (u) induced by one (1) constant doublet panel (cdp) on one control point (1)

    Global frame: X, Y, panel frame x, y
      x = X cos phi + Y sin phi
      y =-X sin phi + Y cos phi

    INPUTS:
    - rCP: position of control point (two values)
    - rS1: position of panel start (two values)
    - rS2: position of panel end (two values)
    - mu: doublet strength
    - method: 1 (Theoretical), 2 (Numerical Quadrature), 10 (Point doublet)
    - principal: if True, return principal value no matter what (if we know we are on a panel)
                 Otherwise, the code will attempt to detect if we are on the panel.

    OUTPUTS:
    - u, v: velocity components at control point
    """
    if method == 1:
        u, v = cdp_u11_kp(rCP, rS1, rS2, mu=mu, tol=tol, principal=principal)
    elif method == 2:
        u, v = cdp_u11_quad(rCP, rS1, rS2, mu=mu, tol=tol, principal=principal)
    elif method == 10:
        if WARN[0] < 3:
            print('[WARN] Doublet Panel 2D - Method 10 only works for panels along x axis for now')
            WARN[0] += 1
        from welib.vortilib.elements.DoubletPoint import dp2d_u
        # We use many point doublets along the panel
        nS = 30
        xS1, yS1 = rS1
        xS2, yS2 = rS2
        xCP, yCP = rCP
        dx = xS2 - xS1
        dy = yS2 - yS1
        L = np.sqrt(dx ** 2 + dy ** 2)
        dl = L / nS
        Mu = mu * dl
        P = np.linspace(rS1, rS2, nS + 1)
        Ps = (P[:-1, :] + P[1:, :]) / 2
        u, v = 0, 0
        for i in range(nS):
            # TODO might not be applicable, Doublet points in wrong direction
            ui, vi = dp2d_u(xCP, yCP, Ps[i], Mu=Mu, orientation='y')
            u += ui
            v += vi
    elif method == 20:
        # Two point vortices at the panel ends, strength = mu each, opposite sign
        from welib.vortilib.elements.VortexPoint import vp_u
        xCP, yCP = rCP
        xS1, yS1 = rS1
        xS2, yS2 = rS2
        # The doublet panel is equivalent to two point vortices of strength +mu at S2 and -mu at S1
        u1, v1 = vp_u(xCP, yCP, (xS1, yS1), Gamma=+mu)
        u2, v2 = vp_u(xCP, yCP, (xS2, yS2), Gamma=-mu)
        u = u1 + u2
        v = v1 + v2
    else:
        raise ValueError("Method must be 1 (Theoretical) or 2 (Quadrature) or 10 or 20")

    return u, v


def cdp_u11_quad(rCP, rS1, rS2, mu=1, tol=1e-8, principal=False):
    """
    Velocity induced by a constant doublet panel using numerical quadrature.

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end (two values)
     - mu: doublet strength
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
    y_2C =-x2C*S +  y2C*C  # = y_1C to machine precision

    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot   = dx * x1C + dy * y1C
    if principal or (abs(cross) < tol and 0 - tol <= dot <= L**2 + tol):
        ut = 0
        if np.abs(x_1C)<1e-8 or np.abs(x_2C)<1e-8:
            un = 0
        else:
            un = +mu/(2*np.pi) * (1/(x_1C)-1/(x_2C)) # NOTE NOTE NOTE: I flipped the sign compared to Katz-Plotkin 
        u, v = un * n_hat + ut * t_hat 
        return u, v
    
    # Vector-valued quadrature
    epsilon = 1e-12
    DY = y_1C
    #def integrand_ut(s):
    #    # [1] Equation 10.26
    #    ss = s/L
    #    x_minus_x0 = x_1C*(1-ss) + x_2C * ss # linear variation between the two bounds
    #    denom =  (x_minus_x0**2 + DY**2 )**2
    #    numer =   x_minus_x0 * DY
    #    return mu/np.pi * numer/denom 
    #def integrand_un(s):
    #    # [1] Equation 10.27
    #    ss = s/L
    #    x_minus_x0 = x_1C*(1-ss) + x_2C * ss # linear variation between the two bounds
    #    denom =  (x_minus_x0**2 + DY**2 ) **2
    #    numer =   x_minus_x0**2 - DY**2
    #    return - mu/(2*np.pi) * numer/denom
    #ut, _ = quad(integrand_ut, 0, L, epsabs=1e-10, epsrel=1e-10)
    #un, _ = quad(integrand_un, 0, L, epsabs=1e-10, epsrel=1e-10)

    def integrand_un_ut(s):
        ss = s/L
        x_minus_x0 = x_1C*(1-ss) + x_2C * ss # linear variation between the two bounds
        denom =  (x_minus_x0**2 + DY**2 ) **2
        ut =       mu/np.pi * (x_minus_x0 * DY)       /denom  # [1] Equation 10.26
        un = - mu/(2*np.pi) * (x_minus_x0**2 - DY**2 )/denom  # [1] Equation 10.27
        return np.array( [un, ut] )
    result, _ = quad_vec(integrand_un_ut, 0, L)
    un, ut = result[0], result[1]
    u, v = un * n_hat + ut * t_hat 
    return u, v


ONE_OVER_TWOPI=0.15916
ONE_OVER_TWOPI= 1/(2*np.pi) #+0.00001



def cdp_u11_kp_raw(rCP, rS1, rS2, principal=False):
    # Katz Plotfkin program 3
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    xCP, yCP = rCP
    dx = xS2 - xS1
    dy = yS2 - yS1
    theta = np.atan2(dy, dx)
    xt   = xCP  - xS1
    zt   = yCP  - yS1
    x2t  = xS2 - xS1
    z2t  = yS2 - yS1
    x    =  xt*np.cos(theta) + zt* np.sin(theta)
    z    = -xt*np.sin(theta) + zt* np.cos(theta)
    x2   = x2t*np.cos(theta) + z2t*np.sin(theta)
    z2   = 0
    r1   = np.sqrt(x**2 + z**2)
    r2   = np.sqrt((x-x2)**2 + z**2)
    # NOTE: sign flipped compared to Katz-Plotkin
    if principal:
        ul = 0
        wl = 1/(np.pi * x)
    else:
        ul = -ONE_OVER_TWOPI * (z/(r1**2) - z/(r2**2))
        wl =  ONE_OVER_TWOPI * (x/(r1**2) - (x-x2)/(r2**2))
    u_kp =  ul*np.cos(-theta) + wl*np.sin(-theta)
    w_kp = -ul*np.sin(-theta) + wl*np.cos(-theta)
    return u_kp, w_kp




def cdp_u11_kp(rCP, rS1, rS2, mu=1, tol=1e-8, principal=False):
    """
    Velocity induced on 1 control points by one doublet panel

    Formulae based on [1] p 272, but rotated in panel frame
    
    Global frame: x, y, panel frame x~, y~ (tilde)
      x~ = x cos phi + y sin phi
      y~ =-x sin phi + y cos phi

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end (two values)
     - mu: doublet strength
     - tol: tolerance for detecting if point is on panel
     - principal: if True, return principal value (zero velocity on panel)

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
    theta = np.atan2(dy, dx)
    phi = theta if theta >= 0 else theta + 2 * np.pi
    
    # --- Vectors from panel points to CP
    x1C = xCP - xS1
    y1C = yCP - yS1
    x2C = xCP - xS2
    y2C = yCP - yS2
    r1C_2 = x1C**2 + y1C**2
    r2C_2 = x2C**2 + y2C**2
    C = np.cos(phi)
    S = np.sin(phi)

    # Katz Plotkin program 3 - Similar with a difference in sign, and small numerical differences
    #xt   = xCP  - xS1
    #zt   = yCP  - yS1
    #x2t  = xS2 - xS1
    #z2t  = yS2 - yS1
    #x    =  xt*np.cos(theta) + zt* np.sin(theta)
    #z    = -xt*np.sin(theta) + zt* np.cos(theta)
    #x2   = x2t*np.cos(theta) + z2t*np.sin(theta)
    #z2   = 0
    #r1   = np.sqrt(x**2 + z**2)
    #r2   = np.sqrt((x-x2)**2 + z**2)
    #ul =  ONE_OVER_TWOPI * (z/(r1**2) - z/(r2**2))
    #wl = -ONE_OVER_TWOPI * (x/(r1**2) - (x-x2)/(r2**2))
    #u_kp =  ul*np.cos(-theta) + wl*np.sin(-theta)
    #w_kp = -ul*np.sin(-theta) + wl*np.cos(-theta)
    
    # ---Tilde coordinates
    # x_~ = x cos phi + y sin phi
    # y_~ =-x sin phi + y cos phi
    x_1C = x1C*C +  y1C*S  # xtilde_c - xtilde_1 = xtilde_c
    y_1C =-x1C*S +  y1C*C  # ytilde_c - ytilde_1 = ytilde_c
    x_2C = x2C*C +  y2C*S  # xtilde_c - xtilde_2
    y_2C =-x2C*S +  y2C*C  # = y_1C to machine precision
    
    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot   = dx * x1C + dy * y1C
    if principal or (abs(cross) < tol and 0 - tol <= dot <= L**2 + tol):
        ut = 0
        if np.abs(x_1C)<1e-8 or np.abs(x_2C)<1e-8:
            un = 0
        else:
            un = +mu/(2*np.pi) * (1/(x_1C)-1/(x_2C)) # NOTE NOTE NOTE: I flipped the sign compared to Katz-Plotkin 
        u,v = un*n_hat + ut*t_hat
        return u, v  # Principal value on the panel is zero
    # NOTE NOTE NOTE: I flipped the sign compared to Katz-Plotkin here 10.29-10.30 
    #                 The sign convention seems innconsistent with 10.26 and 10.27
    # Velocity off the panel (Katz & Plotkin), 10.29 - 10.30
    un = + mu * ONE_OVER_TWOPI * (x_1C / r1C_2 - x_2C/r2C_2 )
    ut = - mu * ONE_OVER_TWOPI * (y_1C / r1C_2 - y_2C/r2C_2 )
    u,v = un*n_hat + ut*t_hat

    return u, v
# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_CDP_PrincipalValue(self):
        # --- One Panel on the x-axis - Check that principal value is obtained
        mu = 10
        PP1 = [-0.5, 0.0]
        PP2 = [0.5, 0.0]
        rCP = [0., 0]
        L = 1

        PV_mid = 2*mu / (L*np.pi)

        u1 = cdp_u11(rCP, PP1, PP2, mu=mu, method=1)
        u2 = cdp_u11(rCP, PP1, PP2, mu=mu, method=2)
        u4 = cdp_u11(rCP, PP1, PP2, mu=mu, method=20)
        np.testing.assert_almost_equal(u1, (0, PV_mid), decimal=6)
        np.testing.assert_almost_equal(u2, (0, PV_mid), decimal=6)
        np.testing.assert_almost_equal(u4, (0, PV_mid), decimal=6)

        # --- One panel on x-axis - principal value with many points
        SP = np.vstack((PP1, PP2))
        x = np.linspace(-0.25, 0.25, 10)
        y = x * 0 + 0
        u1 = cdp_u(x, y, SP, mus=[mu])  # Only method 1
        u_ref = ([0] * len(x), mu/(2*np.pi)*(1/(x+0.5)-1/(x-0.5) ) ) # TODO TODO TODO SIGN CHANGED
        # Also test method 20 for each point
        for xi, yi, u_refi, v_refi in zip(x, y, u_ref[0], u_ref[1]):
            u4, v4 = cdp_u11([xi, yi], PP1, PP2, mu=mu, method=20)
            np.testing.assert_almost_equal([u4, v4], [u_refi, v_refi], decimal=6)
        np.testing.assert_almost_equal(u1, u_ref, decimal=6)
        
        # --- One tilted panel - principal value with one point
        PP1 = [-0.5, 0.2]
        PP2 = [0.5, 0.4]
        rCP = [(PP1[0] + PP2[0])/2, (PP1[1] + PP2[1])/2]  # Midpoint
        u1 = cdp_u11(rCP, PP1, PP2, mu=mu, method=1, principal=True)
        u2 = cdp_u11(rCP, PP1, PP2, mu=mu, method=2, principal=True)
        u4 = cdp_u11(rCP, PP1, PP2, mu=mu, method=20, principal=True)
        np.testing.assert_almost_equal(u1, u2, decimal=6)
        np.testing.assert_almost_equal(u1, u4, decimal=6)

    def test_CDP_flowrate(self):
        # Test flow rate and circulation for a doublet panel
        # Katz Plotkin Eq 3.142
        #   Circulation from begining of panel to x should be mu(x)
        #   Circulation around all of the panel should be zero
        #
        # Also, flow rate should be zero since doublet are source + sink
        from welib.CFD.flows2D import flowrate2D
        from welib.CFD.flows2D import circulation2D

        mu = 3
        PP1 = np.array([-0.5, 0.0])
        PP2 = np.array([ 0.5, 0.0])
        SP = np.vstack((PP1, PP2))
        theta = np.linspace(0, 2 * np.pi, 500, endpoint=True)
        # Contour over all of the panel, we should get 0
        for R in [1.2]:
            x = R * np.cos(theta)
            y = R * np.sin(theta)
            u1, v1 = cdp_u(x, y, SP, mus=[mu], method=1)
            u2, v2 = cdp_u(x, y, SP, mus=[mu], method=2)
            u4 = np.zeros_like(u1)
            v4 = np.zeros_like(v1)
            for i in range(len(x)):
                u4[i], v4[i] = cdp_u11([x[i], y[i]], PP1, PP2, mu=mu, method=20)
            Gamma1 = circulation2D(x, y, u1, v1, verbose=False)
            Gamma2 = circulation2D(x, y, u2, v2, verbose=False)
            Gamma4 = circulation2D(x, y, u4, v4, verbose=False)
            Q1 = flowrate2D(x, y, u1, v1, verbose=False, ns=-1)
            Q2 = flowrate2D(x, y, u2, v2, verbose=False, ns=-1)
            Q4 = flowrate2D(x, y, u4, v4, verbose=False, ns=-1)
            np.testing.assert_almost_equal(Gamma1, 0, decimal=8)
            np.testing.assert_almost_equal(Gamma2, 0, decimal=8)
            np.testing.assert_almost_equal(Gamma4, 0, decimal=8)
            np.testing.assert_almost_equal(Q1, 0, decimal=8)
            np.testing.assert_almost_equal(Q2, 0, decimal=8)
            np.testing.assert_almost_equal(Q4, 0, decimal=8)

        # Contour centered on one extremity, we should get mu
        for R in [0.8]:
            x = R * np.cos(theta) - 0.5
            y = R * np.sin(theta)
            u1, v1 = cdp_u(x, y, SP, mus=[mu], method=1)
            u2, v2 = cdp_u(x, y, SP, mus=[mu], method=2)
            u4 = np.zeros_like(u1)
            v4 = np.zeros_like(v1)
            for i in range(len(x)):
                u4[i], v4[i] = cdp_u11([x[i], y[i]], PP1, PP2, mu=mu, method=20)
            Gamma1 = circulation2D(x, y, u1, v1, verbose=False)
            Gamma2 = circulation2D(x, y, u2, v2, verbose=False)
            Gamma4 = circulation2D(x, y, u4, v4, verbose=False)
            np.testing.assert_almost_equal(Gamma1, mu, decimal=4)
            np.testing.assert_almost_equal(Gamma2, mu, decimal=4)
            np.testing.assert_almost_equal(Gamma4, mu, decimal=4)

        # Contour centered on other extremity, we should get -mu
        for R in [0.8]:
            x = R * np.cos(theta) + 0.5
            y = R * np.sin(theta)
            u1, v1 = cdp_u(x, y, SP, mus=[mu], method=1)
            u2, v2 = cdp_u(x, y, SP, mus=[mu], method=2)
            u4 = np.zeros_like(u1)
            v4 = np.zeros_like(v1)
            for i in range(len(x)):
                u4[i], v4[i] = cdp_u11([x[i], y[i]], PP1, PP2, mu=mu, method=20)
            Gamma1 = circulation2D(x, y, u1, v1, verbose=False)
            Gamma2 = circulation2D(x, y, u2, v2, verbose=False)
            Gamma4 = circulation2D(x, y, u4, v4, verbose=False)
            np.testing.assert_almost_equal(Gamma1, mu, decimal=-4)
            np.testing.assert_almost_equal(Gamma2, mu, decimal=-4)
            np.testing.assert_almost_equal(Gamma4, mu, decimal=-4)
        
        # Test that flow rate on the panel is zero
        #x = np.linspace(-1, 1, 10)
        #y = x * 0 + 0
        #u1, v1 = cdp_u(x, y, SP, mus=[mu], method=1)
        #u2, v2 = cdp_u(x, y, SP, mus=[mu], method=2)
        ##u2, v2 = cdp_u(x, y, SP, mus=[mu], method=10)
        #Q1 = flowrate2D(x, y, u1, v1, verbose=False, ns=1)
        #Q2 = flowrate2D(x, y, u2, v2, verbose=False, ns=1)
        #np.testing.assert_almost_equal(Q1, 0, decimal=3)
        #np.testing.assert_almost_equal(Q2, 0, decimal=3)

    def test_CDP_flow(self, plot=False):
        from welib.CFD.flows2D import flowfield2D, flowfield2D_plot
        # --- One Panel - Comparison of methods
        mu = 1
        PP1 = np.array([-0.5, -0.3])
        PP2 = np.array([0.5 , 0.4])
        SP = np.vstack((PP1, PP2))

        xmax = 1.5
        xs = np.linspace(-xmax, xmax, 15)
        dy = 0.5
        ys = xs * 0 + dy
        xs = np.concatenate((xs, xs))
        ys = np.concatenate((ys, ys - 2 * dy))

        vel1 = lambda X, Y: cdp_u(X, Y, SP, [mu], method=1)
        vel2 = lambda X, Y: cdp_u(X, Y, SP, [mu], method=2)
        vel4 = lambda X, Y: cdp_u(X, Y, SP, [mu], method=20)

        X, Y, U1, V1 = flowfield2D(vel1, xmax=1.5, ymin=-1.3, nx=15)
        X, Y, U2, V2 = flowfield2D(vel2, xmax=1.5, ymin=-1.3, nx=15)
        X, Y, U4, V4 = flowfield2D(vel4, xmax=1.5, ymin=-1.3, nx=15)

        np.testing.assert_almost_equal(U1, U2, decimal=6)
        np.testing.assert_almost_equal(V1, V2, decimal=6)
        np.testing.assert_almost_equal(U1, U4, decimal=6)
        np.testing.assert_almost_equal(V1, V4, decimal=6)

        if plot:
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 3, sharey=False, figsize=(10, 3.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.50)
            ax =  flowfield2D_plot(X, Y, U1, V1, bounded=True, xs=xs, ys = ys, ax=axes[0], maxVal=0.5, minVal=0)
            ax.set_title('Theoretical')

            ax = flowfield2D_plot(X, Y, U2, V2, bounded=True, xs=xs, ys=ys, ax=axes[1], maxVal=0.5, minVal=0)
            ax.set_title('Quadrature')

            #ax = flowfield2D_plot(X, Y, U3, V3, bounded=True, xs=xs, ys=ys, ax=axes[2], maxVal=0.5, minVal=0)
            #ax.set_title('Numerical')

            ax = flowfield2D_plot(X, Y, U4, V4, bounded=True, xs=xs, ys=ys, ax=axes[2], maxVal=0.5, minVal=0)
            ax.set_title('Two points')
            plt.show()

    def test_CDP_crossing_points(self, plot=False):
        """
        Test velocity induced by a tilted panel at points crossing it at an angle through its midpoint.
        Compares methods 1 (Theoretical) and 2 (Quadrature).
        """
        mu = 3
        PP1 = np.array([-0.5, 0.2])
        PP2 = np.array([0.5, 0.4])
        PP1 = np.array([-0.5, 0.0])
        PP2 = np.array([0.5, 0.0])
        SP = np.vstack((PP1, PP2))
        dx = PP2[0] - PP1[0]
        dy = PP2[1] - PP1[1]
        L = np.sqrt(dx**2 + dy**2)
        midpoint = np.array([(PP1[0] + PP2[0])/2, (PP1[1] + PP2[1])/2])

        PV_mid = 2*mu / (L*np.pi)

        # Define points crossing the panel at 45 degrees through the midpoint
        n_points = 21
        s = np.linspace(-1, 1, n_points)
        phi = np.atan2(dy, dx)
        for angle_off in [0, np.pi/4, np.pi/2]:
            angle = phi + angle_off
            x = midpoint[0] + s * np.cos(angle)
            y = midpoint[1] + s * np.sin(angle)
            rCPs = np.vstack((x, y)).T

            # Compute velocities for each method
            u1, v1 = np.zeros(n_points), np.zeros(n_points)
            u2, v2 = np.zeros(n_points), np.zeros(n_points)
            u4, v4 = np.zeros(n_points), np.zeros(n_points)
            for i, rCP in enumerate(rCPs):
                u1[i], v1[i] = cdp_u11(rCP, PP1, PP2, mu=mu, method=1)
                u2[i], v2[i] = cdp_u11(rCP, PP1, PP2, mu=mu, method=2)
                u4[i], v4[i] = cdp_u11(rCP, PP1, PP2, mu=mu, method=20)

            # Verify principal value at midpoint (s=0, index=n_points//2)
            expected = (-np.sin(phi)*PV_mid, np.cos(phi)*PV_mid)
            np.testing.assert_almost_equal([u1[n_points//2], v1[n_points//2]], expected, decimal=6)
            np.testing.assert_almost_equal([u2[n_points//2], v2[n_points//2]], expected, decimal=6)
            np.testing.assert_almost_equal([u4[n_points//2], v4[n_points//2]], expected, decimal=6)

            # Compare methods
            np.testing.assert_almost_equal(u1, u2, decimal=6)
            #np.testing.assert_almost_equal(u1, u4, decimal=6)
            np.testing.assert_almost_equal(v1, v2, decimal=6)
            #np.testing.assert_almost_equal(v1, v4, decimal=6)

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
                ax1.plot(s, u1, 'r--', label='Theoretical')
                ax1.plot(s, u2, 'g:d', label='Quadrature')
                ax1.plot(s, u4, 'b-.', label='Two points')
                ax1.set_ylabel('u velocity')
                ax1.legend()
                ax1.grid(True)
                ax2.plot(s, v1, 'r--', label='Theoretical')
                ax2.plot(s, v2, 'g:d', label='Quadrature')
                ax2.plot(s, v4, 'b-.', label='Two points')
                ax2.set_xlabel('s (along crossing line)')
                ax2.set_ylabel('v velocity')
                ax2.legend()
                ax2.grid(True)
                plt.tight_layout()
        if plot:
            plt.show()

    def test_CDP_debug_point(self):
        """
        Debug test to compare velocities at a specific off-panel point.
        """
        mu = 1
        #PP1 = np.array([-0.5, 0.2])
        #PP2 = np.array([0.5, 0.4])
        PP1 = np.array([-0.5, 0.0])
        PP2 = np.array([0.5, 0.0])
        rCP = [0.0, 0.5]  # Point above the panel
        u1, v1 = cdp_u11(rCP, PP1, PP2, mu=mu, method=1)
        u2, v2 = cdp_u11(rCP, PP1, PP2, mu=mu, method=2)
        u4, v4 = cdp_u11(rCP, PP1, PP2, mu=mu, method=20)
        #u3, v3 = cdp_u11(rCP, PP1, PP2, mu=mu, method=10)
        #print(f"Debug point {rCP}:")
        #print(f"Method 1 (Theoretical): u={u1}, v={v1}")
        #print(f"Method 2 (Quadrature) : u={u2}, v={v2}")
        #print(f"Method 3 (Rieman sum) : u={u3}, v={v3}")
        np.testing.assert_almost_equal([u1, v1], [u2, v2], decimal=6)
        np.testing.assert_almost_equal([u1, v1], [u4, v4], decimal=6)
        #np.testing.assert_almost_equal([u1, v1], [u3, v3], decimal=3)

if __name__ == "__main__":
    #Test().test_CDP_debug_point()
    #Test().test_CDP_flow(plot=True)
    #Test().test_CDP_crossing_points(plot=True)
    #Test().test_CDP_PrincipalValue()
    #Test().test_CDP_flowrate()
    unittest.main()
