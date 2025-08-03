""" 
Constant Source Panel 2D (CSP)

Reference: 
 [0] Anderson Fundamentals-of-aerodynamics-6-Edition, p.292
 [1] Katz - Plotkin - Low speed aerodynamics p 270
"""
import numpy as np
import unittest
from scipy.integrate import quad, quad_vec
# --------------------------------------------------------------------------------}
# --- Wrappers 
# --------------------------------------------------------------------------------{
def dcsp_u(X, Y, SP1, SP2, sigmas, debug=False):
    """ 
    Constant Source Panels, delimited by point 1 and 2 (can be discontinuous)

    INPUTS:
     - X: nd-array of control point x-coordinates
     - Y: nd-array of control point y-coordinates
     - SP1: (n,2)-array of first source panels points
     - SP2: (n,2)-array of second source panels points
     - sigmas: (n)-array of source panels intensities
    """
    X = np.asarray(X)
    Y = np.asarray(Y)
    shp = X.shape
    X = X.flatten()
    Y = Y.flatten()
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    for j in range(len(SP1)): # Loop on panels
        if debug:
            print('Panel',j, SP1[j], SP2[j], sigmas[j])
        ui, vi = csp_u1N(X, Y, SP1[j], SP2[j], sigmas[j], debug=debug)
        U+=ui
        V+=vi
    U = U.reshape(shp)
    V = V.reshape(shp)
    return U, V


def csp_u(X, Y, SP, sigmas, method=1, debug=False):
    """ 
    Contiguous Constant Source Panels 
    Continguous => Panels are formed by consecutive points, can potentially form a closed loop

    INPUTS:
     - X: nd-array of control point x-coordinates
     - Y: nd-array of control point y-coordinates
     - SP: (n,2)-array of source panels points
     - sigmas: (n-1)-array of source panels intensities
    """
    X = np.asarray(X)
    Y = np.asarray(Y)
    shp = X.shape
    X = X.flatten()
    Y = Y.flatten()
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    if method==1:
        # Optimized for many control points
        for j in range(len(SP)-1): # Loop on panels
            ui, vi = csp_u1N(X, Y, SP[j], SP[j+1], sigmas[j], debug=debug)
            U+=ui
            V+=vi
    else:
        # Not optimized, for any method
        for i in range(len(X)):
            for j in range(len(SP)-1): # Loop on panels
                ui, vi = csp_u11((X[i],Y[i]), SP[j], SP[j+1], sigmas[j], method=method)
                U[i]+=ui
                V[i]+=vi
    U = U.reshape(shp)
    V = V.reshape(shp)
    return U, V

# --------------------------------------------------------------------------------}
# --- Vectorized Function for One Panel, Many Points
# --------------------------------------------------------------------------------{
def csp_u1N(xCP, yCP, rS1, rS2, sigma=1, tol=1e-8, debug=False):
    """
    Velocity induced on N control points by one constant source panel (csp)

    Formulae based on [1] p 268, but rotated in panel frame

    Same as csp_u11_kp ("method 1", Katz-Plotkin) but vectorized.

    INPUTS:
     - xCP, yCP : position of control points , n-array
     - rS1: position of panel start (2 values)
     - rS2: position of panel end   (2 values)
    OUTPUTS:
     - U, V: velocities at control points, n-array
    """
    U = np.zeros_like(xCP)
    V = np.zeros_like(yCP)
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    
    # --- Panel
    dx = xS2 - xS1
    dy = yS2 - yS1
    L = np.sqrt(dx**2 + dy**2)
    n_hat = np.array([-dy/L, dx/L])
    t_hat = np.array([dx/L , dy/L])
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
    # x_ = X cos phi + Y sin phi
    # y_ =-X sin phi + Y cos phi
    x_1C = x1C*C +  y1C*S 
    y_1C =-x1C*S +  y1C*C
    x_2C = x2C*C +  y2C*S
    y_2C =-x2C*S +  y2C*C
    
    # --- Angles measured in panel frame
    theta1 = np.arctan2(y_1C, x_1C)
    theta2 = np.arctan2(y_2C, x_2C)
    theta1[theta1<0] +=2*np.pi
    theta2[theta2<0] +=2*np.pi
    dtheta = theta2 - theta1

    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot = dx * x1C + dy * y1C
    bOnLine = abs(cross)<tol
    bOnSeg = np.logical_and(- tol <= dot, dot <= L**2 + tol)
    b = np.logical_and(bOnLine, bOnSeg)
    U[b] =  0.5 * sigma * n_hat[0]
    V[b] =  0.5 * sigma * n_hat[1]
    Un = sigma / (2 * np.pi) * dtheta[~b]
    Ut = sigma / (4 * np.pi) * np.log(r1C_2[~b]/r2C_2[~b])
    U[~b] = Un*n_hat[0] + Ut*t_hat[0]
    V[~b] = Un*n_hat[1] + Ut*t_hat[1]
    if debug:
        print('xCP ', xCP, yCP)
        print('xS1 ', xS1, yS1)
        print('x1C', x1C, y1C)
        print('L, dx, dy', L, dx, dy)
        print('cross, dot', cross, dot)
        print('n_hat', n_hat, sigma, b)
    return U, V

# --------------------------------------------------------------------------------}
# --- 11 : Velocity (u) induced by one (1) constant source panel (csp) on one control point (1)
# --------------------------------------------------------------------------------{
def csp_u11(rCP, rS1, rS2, sigma=1, method=1, tol=1e-8, principal=False):
    """
    Velocity (u) induced by one (1) constant source panel (csp) on one control point (1)

    Global frame: X, Y, panel frame x, y
      x = X cos phi + Y sin phi
      y =-X sin phi + Y cos phi

    INPUTS:
    - rCP: position of control point (two values)
    - rS1: position of panel start (two values)
    - rS2: position of panel end (two values)
    - sigma: source strength
    - method: 0 (Anderson), 1 (Katz & Plotkin), 2 (Numerical Quadrature), 10 (Point sources)
    - principal: if True, return principal value no matter what (if we know we are on a panel)
                 Otherwise, the code will attempt to detect if we are on the panel.

    OUTPUTS:
    - u, v: velocity components at control point                 
    """
    if method==0:
        u,v = csp_u11_anderson(rCP, rS1, rS2, sigma=sigma, tol=tol, principal=principal)
    elif method ==1:
        u, v = csp_u11_kp(rCP, rS1, rS2, sigma=sigma, tol=tol, principal=principal)
    elif method == 2:
        u, v = csp_u11_quad(rCP, rS1, rS2, sigma=sigma, tol=tol, principal=principal)
    elif method ==10:
        from welib.vortilib.elements.SourcePoint import sp2d_u
        # We use point many point sources along the panel
        nS = 30
        xS1, yS1 = rS1
        xS2, yS2 = rS2
        xCP, yCP = rCP
        dx = xS2 - xS1
        dy = yS2 - yS1
        L = np.sqrt(dx**2 + dy**2)
        dl = L/nS
        Sig = sigma*dl
        P = np.linspace(rS1, rS2, nS+1)
        Ps  = (P[:-1,:] + P[1:,:]) / 2
        u,v=0,0
        for i in range(nS):
            ui, vi = sp2d_u(xCP, yCP, Ps[i], Sigma=Sig)
            u+=ui
            v+=vi
    else:
        raise ValueError("Method must be 0 (Anderson), 1 (Theoretical), 2 (Quadrature) or 10 (Points)")
            
    return u, v


def csp_u11_quad(rCP, rS1, rS2, sigma=1, tol=1e-8, principal=False):
    """
    Velocity induced by a constant source panel at a control point using numerical quadrature.

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end (two values)
     - sigma: source strength
     - tol: tolerance for detecting if point is on panel
     - principal: if True, return principal value (normal velocity = sigma/2)

    OUTPUTS:
     - u, v: velocity components at control point
    """
    xCP, yCP = rCP
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    dx = xS2 - xS1
    dy = yS2 - yS1
    length = np.sqrt(dx**2 + dy**2)
    n_hat = np.array([-dy/length, dx/length])  # Normal vector

    # Check if point is on the panel
    x1C = xCP - xS1
    y1C = yCP - yS1
    cross = dx * y1C - dy * x1C
    dot = dx * x1C + dy * y1C
    if principal or (abs(cross) < tol and 0 - tol <= dot <= length**2 + tol):
        return 0.5 * sigma * n_hat  # Principal value

    # Numerical quadrature
    epsilon = 1e-10  # Prevent division by zero
    # KEEP ME
    #def integrand(s, component):
    #    xp = xS1 + s * dx / length
    #    yp = yS1 + s * dy / length
    #    r = np.sqrt((xCP - xp)**2 + (yCP - yp)**2 + epsilon)
    #    if component == 'u':
    #        return (sigma / (2 * np.pi)) * (xCP - xp) / r**2
    #    else:
    #        return (sigma / (2 * np.pi)) * (yCP - yp) / r**2
    #u, _ = quad(lambda s: integrand(s, 'u'), 0, length)
    #v, _ = quad(lambda s: integrand(s, 'v'), 0, length)

    def integrand(s):
        xp = xS1 + s * dx / length
        yp = yS1 + s * dy / length
        r = np.sqrt((xCP - xp)**2 + (yCP - yp)**2 + epsilon)
        factor = sigma / (2 * np.pi * r**2)
        return np.array([factor * (xCP - xp), factor * (yCP - yp)])

    result, _ = quad_vec(integrand, 0, length)

    return result[0], result[1]

def csp_u11_kp(rCP, rS1, rS2, sigma=1, tol=1e-8, principal=False):
    """
    Velocity induced on 1 control points by one cource panel

    Formulae based on [1] p 268, but rotated in panel frame

    Global frame: x, y, panel frame x~, y~ (tilde)
      x~ = x cos phi + y sin phi
      y~ =-x sin phi + y cos phi

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end   (two values)
     - sigma: source strength
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
    t_hat = np.array([dx/L , dy/L])
    if principal:
        return 0.5 * sigma * n_hat # On the segment, return principal value
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
    x_1C = x1C*C +  y1C*S  
    y_1C =-x1C*S +  y1C*C
    x_2C = x2C*C +  y2C*S
    y_2C =-x2C*S +  y2C*C  # = y_1C to machine precision
       
    # --- Angles measured in panel frame
    theta1 = np.arctan2(y_1C, x_1C)
    theta2 = np.arctan2(y_2C, x_2C)
    if theta2<0: 
        theta2 += 2*np.pi # Panel angles are assumed positive
    if theta1<0: 
        theta1 += 2*np.pi # Panel angles are assumed positive
    dtheta = theta2 - theta1

    # --- Check if the point is on the line using cross product
    cross = dx * y1C - dy * x1C
    dot   = dx * x1C + dy * y1C
    if abs(cross) < tol: # We are on the line or within the segment, return PV
        if 0 - tol <= dot <= L**2 + tol:
            u,v = 0.5 * sigma * n_hat # On the panel, return principal value
            return u,v
    un = sigma / (2 * np.pi) * dtheta 
    ut = sigma / (4 * np.pi) * np.log(r1C_2/r2C_2)
    u,v = un*n_hat + ut*t_hat
    return u,v


def csp_u11_anderson(rCP, rS1, rS2, sigma=1, tol=1e-8, principal=False):
    """
    Velocity induced on one control point by one panel
    Reference [1]

    INPUTS:
     - rCP: position of control point (two values)
     - rS1: position of panel start (two values)
     - rS2: position of panel end   (two values)
    OUTPUTS:
     - u,v : velocity at control point
    """
    xS1, yS1 = rS1
    xS2, yS2 = rS2
    xCP, yCP = rCP
    # Panel geometry
    dx = xS2 - xS1
    dy = yS2 - yS1
    ds = np.sqrt(dx**2 + dy**2)
    phi = np.atan2(dy, dx)
    phi = phi if phi>=0 else phi+2*np.pi # positive angles # phi=mod(phi, 2*pi)
    # --- Anderson formulation
    x1C = xCP-xS1
    y1C = yCP-yS1
    CX = -np.cos(phi)
    CY = -np.sin(phi)
    if principal:
        return  CY * 0.5 * sigma, -CX * 0.5 * sigma # -sin phi = nx ,  cos phi = ny

    A  = x1C*CX + y1C*CY
    B  = x1C**2 + y1C**2
    E2 = B-A**2
    # --- Check if point on the line/segment
    cross = dx * y1C - dy * x1C # Cross product
    dot   = dx * x1C + dy * y1C
    if abs(cross) < tol: # We are on the line or within the segment, return PV
        if 0 - tol <= dot <= ds**2 + tol:
            return  CY * 0.5 * sigma, -CX * 0.5 * sigma # -sin phi = nx ,  cos phi = ny
        else:
            print('>>> On line', xCP, yCP) # Mystery case
            return 0,0
    # --- 
    E  = np.sqrt(E2)
    # See [1] Eq. 3.163 adapted for x and y
    LOG   = np.log((ds**2 + 2*A*ds+B)/B) * 0.5
    ANGLE = np.atan2((ds+A),E)-np.atan2(A,E)
    u     = CX*LOG + ((x1C-A*CX)/E)*ANGLE
    v     = CY*LOG + ((y1C-A*CY)/E)*ANGLE
    u *= sigma/(2*np.pi)
    v *= sigma/(2*np.pi)
    return u, v

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_CSP_PrincipalValue(self):
        # --- One Panel on the x axis - Check that principal value is obtained
        sigma = 10
        PP1   = [-0.5,0.0]
        PP2   = [ 0.5,0.0]
        rCP = [0.,0]
        u0 = csp_u11(rCP, PP1, PP2, sigma=sigma, method=0)
        u1 = csp_u11(rCP, PP1, PP2, sigma=sigma, method=1)
        u2 = csp_u11(rCP, PP1, PP2, sigma=sigma, method=2)
        np.testing.assert_almost_equal(u0, (0, 5))
        np.testing.assert_almost_equal(u1, (0, 5))
        np.testing.assert_almost_equal(u2, (0, 5))

        # --- One panel on x-axis - principal value with many points
        SP = np.vstack((PP1,PP2))
        x = np.linspace(-0.5, 0.5, 10)
        y = x*0+0
        u0 = csp_u(x, y, SP, sigmas=[sigma]) # Only method 1
        u_ref = ([0]*len(x), [sigma/2]*len(x))
        np.testing.assert_almost_equal(u0, u_ref)
        
        # --- One tilted panel - principal value with one point
        PP1 = [-0.5, 0.2]
        PP2 = [0.5, 0.4]
        dx = PP2[0] - PP1[0]
        dy = PP2[1] - PP1[1]
        L = np.sqrt(dx**2 + dy**2)
        n_hat = np.array([-dy/L, dx/L])
        rCP = [(PP1[0] + PP2[0])/2, (PP1[1] + PP2[1])/2]  # Midpoint
        u0 = csp_u11(rCP, PP1, PP2, sigma=sigma, method=0, principal=True)
        u1 = csp_u11(rCP, PP1, PP2, sigma=sigma, method=1, principal=True)
        u2 = csp_u11(rCP, PP1, PP2, sigma=sigma, method=2, principal=True)
        expected = 0.5 * sigma * n_hat
        np.testing.assert_almost_equal(u0, expected, decimal=6)
        np.testing.assert_almost_equal(u1, expected, decimal=6)
        np.testing.assert_almost_equal(u2, expected, decimal=6)

    def test_CSP_flowrate(self):

        # Test that flow rate for a source panel is sigma*l
        from welib.CFD.flows2D import flowrate2D
        sigma = 3
        PP1   = np.array([-1,0.0])
        PP2   = np.array([ 1,0.0])
        SP = np.vstack((PP1,PP2))
        l = np.linalg.norm(PP1 - PP2)
        theta = np.linspace(0, 2*np.pi, 500, endpoint=True)
        for R in [1.2]:
            x = R* np.cos(theta)
            y = R* np.sin(theta)
            u0, v0 = csp_u(x, y, SP, sigmas=[sigma], method=0)
            u1, v1 = csp_u(x, y, SP, sigmas=[sigma], method=1)
            u2, v2 = csp_u(x, y, SP, sigmas=[sigma], method=2)
            Q0 = flowrate2D(x, y, u0, v0, verbose=False, ns=-1)
            Q1 = flowrate2D(x, y, u1, v1, verbose=False, ns=-1)
            Q2 = flowrate2D(x, y, u2, v2, verbose=False, ns=-1)
            np.testing.assert_almost_equal(Q0, sigma * l, decimal=1)
            np.testing.assert_almost_equal(Q1, sigma * l, decimal=3)
            np.testing.assert_almost_equal(Q2, sigma * l, decimal=3)
        
        # Test that flow rate on the panel is half sigma l
        x= np.linspace(-1, 1, 10)
        y= x*0+0
        u0, v0 = csp_u(x, y, SP, sigmas=[sigma], method=0)
        u1, v1 = csp_u(x, y, SP, sigmas=[sigma], method=1)
        u2, v2 = csp_u(x, y, SP, sigmas=[sigma], method=2)        
        Q0 = flowrate2D(x, y, u0, v0, verbose=False, ns=1)
        Q1 = flowrate2D(x, y, u1, v1, verbose=False, ns=1)
        Q2 = flowrate2D(x, y, u2, v2, verbose=False, ns=1)
        np.testing.assert_almost_equal(Q0, sigma * l / 2, decimal=3)
        np.testing.assert_almost_equal(Q1, sigma * l / 2, decimal=3)
        np.testing.assert_almost_equal(Q2, sigma * l / 2, decimal=3)

    def test_CSP_flow(self, plot=False):
        from welib.CFD.flows2D import flowfield2D, flowfield2D_plot
        # --- One Panel - Comparison of methods including point sources
        sigma = 1
        PP1 = np.array([-0.5,0])
        PP2 = np.array([ 0.6,0.2])
        SP  = np.vstack((PP1,PP2))

        xmax = 1.5
        xs   = np.linspace(-xmax, xmax, 15)
        dy   = 0.5
        ys   = xs*0 + dy
        xs = np.concatenate((xs,xs))
        ys = np.concatenate((ys,ys-2*dy))

        # --- 
        vel0 = lambda X, Y : csp_u(X, Y, SP, [sigma], method=0)
        vel1 = lambda X, Y : csp_u(X, Y, SP, [sigma], method=1)
        vel2 = lambda X, Y : csp_u(X, Y, SP, [sigma], method=2)
        vel3 = lambda X, Y : csp_u(X, Y, SP, [sigma], method=10)

        X, Y, U0, V0 =  flowfield2D(vel0, xmax=1.5, ymin=0.3, nx=15)
        X, Y, U1, V1 =  flowfield2D(vel1, xmax=1.5, ymin=0.3, nx=15)
        X, Y, U2, V2  = flowfield2D(vel2, xmax=1.5, ymin=0.3, nx=15)
        X, Y, U3, V3 =  flowfield2D(vel3, xmax=1.5, ymin=0.3, nx=15)

        np.testing.assert_almost_equal(U0, U1)
        np.testing.assert_almost_equal(V0, V1)
        np.testing.assert_almost_equal(U0, U2, decimal=6)
        np.testing.assert_almost_equal(V0, V2, decimal=6)        
        np.testing.assert_almost_equal(U0, U3, decimal=3)
        np.testing.assert_almost_equal(V0, V3, decimal=3)

        if plot:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots(1, 4, sharey=False, figsize=(15.4,3.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.50)
            ax =  flowfield2D_plot(X, Y, U0, V0, bounded=False, xs=xs, ys = ys, ax=axes[0], maxVal=0.5, minVal=0)
            ax.set_title('Anderson')
            ax =  flowfield2D_plot(X, Y, U1, V1, bounded=False, xs=xs, ys = ys, ax=axes[1], maxVal=0.5, minVal=0)
            ax.set_title('Katz & Plotkin')
            ax = flowfield2D_plot(X, Y, U2, V2, bounded=False, xs=xs, ys=ys, ax=axes[2], maxVal=0.5, minVal=0)
            ax.set_title('Quadrature')
            ax = flowfield2D_plot(X, Y, U3, V3, bounded=False, xs=xs, ys=ys, ax=axes[3], maxVal=0.5, minVal=0)
            ax.set_title('Point Sources')
            plt.show()


    def test_CSP_crossing_points(self, plot=False):
        """
        Test velocity induced by a tilted panel at points crossing it at an angle through its midpoint.
        Compares methods 0 (Anderson), 1 (Katz & Plotkin), and 2 (Quadrature).
        """
        sigma = 1
        PP1 = np.array([-0.5, 0.2])
        PP2 = np.array([0.5, 0.4])
        SP = np.vstack((PP1, PP2))
        dx = PP2[0] - PP1[0]
        dy = PP2[1] - PP1[1]
        L = np.sqrt(dx**2 + dy**2)
        n_hat = np.array([-dy/L, dx/L])  # Normal vector
        midpoint = np.array([(PP1[0] + PP2[0])/2, (PP1[1] + PP2[1])/2])

        # Define points crossing the panel at 45 degrees through the midpoint
        n_points = 21
        s = np.linspace(-1, 1, n_points)  # Parameter along the crossing line
        # Line direction perpendicular to panel tangent, rotated 45 degrees
        phi = np.atan2(dy, dx)
        angle = phi + np.pi/4  # 45 degrees relative to panel
        x = midpoint[0] + s * np.cos(angle)
        y = midpoint[1] + s * np.sin(angle)
        rCPs = np.vstack((x, y)).T

        # Compute velocities for each method
        u0, v0 = np.zeros(n_points), np.zeros(n_points)
        u1, v1 = np.zeros(n_points), np.zeros(n_points)
        u2, v2 = np.zeros(n_points), np.zeros(n_points)
        for i, rCP in enumerate(rCPs):
            u0[i], v0[i] = csp_u11(rCP, PP1, PP2, sigma=sigma, method=0)
            u1[i], v1[i] = csp_u11(rCP, PP1, PP2, sigma=sigma, method=1)
            u2[i], v2[i] = csp_u11(rCP, PP1, PP2, sigma=sigma, method=2)

        # Verify principal value at midpoint (s=0, index=n_points//2)
        expected = 0.5 * sigma * n_hat
        np.testing.assert_almost_equal([u0[n_points//2], v0[n_points//2]], expected, decimal=6)
        np.testing.assert_almost_equal([u1[n_points//2], v1[n_points//2]], expected, decimal=6)
        np.testing.assert_almost_equal([u2[n_points//2], v2[n_points//2]], expected, decimal=6)

        # Compare methods
        np.testing.assert_almost_equal(u0, u1, decimal=6)
        np.testing.assert_almost_equal(v0, v1, decimal=6)
        np.testing.assert_almost_equal(u0, u2, decimal=6)
        np.testing.assert_almost_equal(v0, v2, decimal=6)

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
            plt.show()

            # Plot 2: u and v velocities
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
            ax1.plot(s, u0, 'b-', label='Anderson')
            ax1.plot(s, u1, 'r--', label='Katz & Plotkin')
            ax1.plot(s, u2, 'g:', label='Quadrature')
            ax1.set_ylabel('u velocity')
            ax1.legend()
            ax1.grid(True)
            ax2.plot(s, v0, 'b-', label='Anderson')
            ax2.plot(s, v1, 'r--', label='Katz & Plotkin')
            ax2.plot(s, v2, 'g:', label='Quadrature')
            ax2.set_xlabel('s (along crossing line)')
            ax2.set_ylabel('v velocity')
            ax2.legend()
            ax2.grid(True)
            plt.tight_layout()
            plt.show()

if __name__ == "__main__":
#     Test().test_CSP_flowrate()
    unittest.main()

