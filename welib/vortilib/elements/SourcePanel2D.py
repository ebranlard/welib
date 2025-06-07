""" 
Constant Source Panel 2D (CSP)

Reference: 
 [1] Anderson Fundamentals-of-aerodynamics-6-Edition, p.292
 [2] Katz - Plotkin - Low speed aerodynamics p 270
"""
import numpy as np
import unittest

# --------------------------------------------------------------------------------}
# --- Wrappers 
# --------------------------------------------------------------------------------{
def csp_u(X, Y, SP1, SP2, sigmas, debug=False):
    """ 
    Constant Source Panels, delimited by point 1 and 2 (can be discontinuous

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


def ccsp_u(X, Y, SP, sigmas, method=1, debug=False):
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


def csp_u1N(xCP, yCP, rS1, rS2, sigma=1, tol=1e-8, debug=False):
    """
    Velocity induced on N control points by one panel

    Formulate based on [2] p 268, but rotated in panel frame

    INPUTS:
     - xCP, yCP : position of control points , n-array
     - rS1: position of panel start (2 values)
     - rS2: position of panel ends  (2 values)
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
    # --- From panel points to CP
    dX_r1 = xCP - xS1
    dY_r1 = yCP - yS1
    dX_r2 = xCP - xS2
    dY_r2 = yCP - yS2
    r12 = dX_r1**2 + dY_r1**2
    r22 = dX_r2**2 + dY_r2**2
    C = np.cos(phi)
    S = np.sin(phi)
    # --- Angles measured in panel frame
    #  x = X cos phi + Y sin phi
    #  y =-X sin phi + Y cos phi
    dx_r1 = dX_r1*C +  dY_r1*S
    dy_r1 =-dX_r1*S +  dY_r1*C
    dx_r2 = dX_r2*C +  dY_r2*S
    dy_r2 =-dX_r2*S +  dY_r2*C
    theta1 = np.arctan2(dy_r1, dx_r1)
    theta2 = np.arctan2(dy_r2, dx_r2)
    theta1[theta1<0] +=2*np.pi
    theta2[theta2<0] +=2*np.pi
    dtheta = theta2 - theta1
    # --- Check if the point is on the line using cross product
    dX_r1 = xCP - xS1
    dY_r1 = yCP - yS1
    cross = dx * dY_r1 - dy * dX_r1 # Cross product 
    dot = dx * dX_r1 + dy * dY_r1
    bOnLine = abs(cross)<tol
    bOnSeg = np.logical_and(- tol <= dot, dot <= L**2 + tol)
    b = np.logical_and(bOnLine, bOnSeg)
    U[b] =  0.5 * sigma * n_hat[0]
    V[b] =  0.5 * sigma * n_hat[1]
    Un = sigma / (2 * np.pi) * dtheta[~b]
    Ut = sigma / (4 * np.pi) * np.log(r12[~b]/r22[~b])
    U[~b] = Un*n_hat[0] + Ut*t_hat[0]
    V[~b] = Un*n_hat[1] + Ut*t_hat[1]
    if debug:
        print('xCP ', xCP, yCP)
        print('xS1 ', xS1, yS1)
        print('dX_r1', dX_r1, dY_r1)
        print('L, dx, dy', L, dx, dy)
        print('cross, dot', cross, dot)
        print('n_hat', n_hat, sigma, b)
    return U, V


# --------------------------------------------------------------------------------}
# --- 11 
# --------------------------------------------------------------------------------{
def csp_u11(rCP, rS1, rS2, sigma=1, method=1, tol=1e-8):
    """
    Global frame: X, Y, panel frame x, y
      x = X cos phi + Y sin phi
      y =-X sin phi + Y cos phi

    INPUTS:
    - rS1: position of panel start (2 values)
    - rS2: position of panel end   (2 values)
    - r : position of control point 
    """
    if method==0:
        u,v = csp_u11_anderson(rCP, rS1, rS2, sigma=sigma, tol=tol)
        
    elif method ==1:
        u, v = csp_u11_kp(rCP, rS1, rS2, sigma=sigma, tol=tol)

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
    return u, v


def csp_u11_kp(rCP, rS1, rS2, sigma=1, tol=1e-8):
    """
    Velocity induced on 1 control points by one panel

    Formulae based on [2] p 268, but rotated in panel frame

    Global frame: X, Y, panel frame x, y
      x = X cos phi + Y sin phi
      y =-X sin phi + Y cos phi

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
    
    # --- Panel
    dx = xS2 - xS1
    dy = yS2 - yS1
    L = np.sqrt(dx**2 + dy**2)
    n_hat = np.array([-dy/L, dx/L])
    t_hat = np.array([dx/L , dy/L])
    phi = np.atan2(dy, dx)
    phi = phi if phi>=0 else phi+2*np.pi # positive angles # phi=mod(phi, 2*pi)

    # See [2] p 268 Rotated in Panel frame
    # --- From panel points to CP
    dX_r1 = xCP - xS1
    dY_r1 = yCP - yS1
    dX_r2 = xCP - xS2
    dY_r2 = yCP - yS2
    r12 = dX_r1**2 + dY_r1**2
    r22 = dX_r2**2 + dY_r2**2
    C = np.cos(phi)
    S = np.sin(phi)
    # --- Angles measured in panel frame
    dx_r1 = dX_r1*C +  dY_r1*S
    dy_r1 =-dX_r1*S +  dY_r1*C
    dx_r2 = dX_r2*C +  dY_r2*S
    dy_r2 =-dX_r2*S +  dY_r2*C
    theta1 = np.arctan2(dy_r1, dx_r1)
    theta2 = np.arctan2(dy_r2, dx_r2)
    if theta2<0: 
        theta2 += 2*np.pi # Panel angles are assumed positive
    if theta1<0: 
        theta1 += 2*np.pi # Panel angles are assumed positive
    dtheta = theta2 - theta1

    # --- Check if the point is on the line using cross product
    dX_r1 = xCP - xS1
    dY_r1 = yCP - yS1
    cross = dx * dY_r1 - dy * dX_r1 # Cross product 
    dot   = dx * dX_r1 + dy * dY_r1
    #print('xCP ', xCP, yCP)
    #print('xS1 ', xS1, yS1)
    #print('dX_r1', dX_r1, dY_r1)
    #print('L, dx, dy', L, dx, dy)
    #print('cross, dot', cross, dot)
    if abs(cross) < tol:
        # We are on the line or within the segment.
        if 0 - tol <= dot <= L**2 + tol:
            #print('>>> On segment ', n_hat, sigma)
            u,v = 0.5 * sigma * n_hat # On the segment, return principal value
            return u,v
    #print('Not on segment')
    un = sigma / (2 * np.pi) * dtheta 
    ut = sigma / (4 * np.pi) * np.log(r12/r22)
    u,v = un*n_hat + ut*t_hat
    return u,v


def csp_u11_anderson(rCP, rS1, rS2, sigma=1, tol=1e-8):
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
    DX1 = xCP-xS1
    DY1 = yCP-yS1
    CX = -np.cos(phi)
    CY = -np.sin(phi)
    A  = DX1*CX + DY1*CY
    B  = DX1**2 + DY1**2
    E2 = B-A**2
    # --- Check if point on the line/segment
    dX_r1 = xCP - xS1
    dY_r1 = yCP - yS1
    cross = dx * dY_r1 - dy * dX_r1 # Cross product
    dot   = dx * dX_r1 + dy * dY_r1
    if abs(cross) < tol:
        # We are on the line or within the segment.
        if 0 - tol <= dot <= ds**2 + tol:
            #print('>>> On segment', xCP, yCP)
            # On the segment, return principal value
            u = CY * 0.5 * sigma # -sin phi  = nx
            v =-CX * 0.5 * sigma #  cos phi  = ny
            return u,v
        else:
            print('>>> On line', xCP, yCP)
            return 0,0
    # --- 
    E  = np.sqrt(E2)
    # See [1] Eq. 3.163 adapted for x and y
    LOG   = np.log((ds**2 + 2*A*ds+B)/B) * 0.5
    ANGLE = np.atan2((ds+A),E)-np.atan2(A,E)
    u     = CX*LOG + ((DX1-A*CX)/E)*ANGLE
    v     = CY*LOG + ((DY1-A*CY)/E)*ANGLE
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
        np.testing.assert_almost_equal(u0, (0, 5))
        np.testing.assert_almost_equal(u1, (0, 5))

        # --- One panel on the x axis, many points on the segment
        SP = np.vstack((PP1,PP2))
        x = np.linspace(-0.5, 0.5, 10)
        y = x*0+0
        u, v = ccsp_u(x, y, SP, sigmas=[sigma], method=1)
        np.testing.assert_almost_equal(u, [0]*len(x))
        np.testing.assert_almost_equal(v, [sigma/2]*len(x))

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
            u, v = ccsp_u(x, y, SP, sigmas=[sigma], method=1)
            #print('x', x)
            #print('y', y)
            #print('u', u)
            #print('v', v)
            Q = flowrate2D(x, y, u, v, verbose=False, ns=-1)
            #print('l', l)
            #print('Q', Q)
            np.testing.assert_almost_equal(Q, sigma*l, 3)
        
        # Test that flow rate on the panel is half sigma l
        x= np.linspace(-1, 1, 10)
        y= x*0+0
        u, v = ccsp_u(x, y, SP, sigmas=[sigma], method=1)
        #print('x', x)
        #print('y', y)
        #print('u', u)
        #print('v', v)
        Q = flowrate2D(x, y, u, v, verbose=False, ns=1)
        #print('l', l)
        #print('Q', Q)
        np.testing.assert_almost_equal(Q, sigma*l/2, 3)

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
#         xs, ys = ys, xs

        # --- 
        vel0 = lambda X, Y : ccsp_u(X, Y, SP, [sigma], method=0)
        vel1 = lambda X, Y : ccsp_u(X, Y, SP, [sigma], method=1)
        vel3 = lambda X, Y : ccsp_u(X, Y, SP, [sigma], method=10)

        X, Y, U0, V0 =  flowfield2D(vel0, xmax=1.5, ymin=0.3, nx=15)
        X, Y, U1, V1 =  flowfield2D(vel1, xmax=1.5, ymin=0.3, nx=15)
        X, Y, U3, V3 =  flowfield2D(vel3, xmax=1.5, ymin=0.3, nx=15)

        np.testing.assert_almost_equal(U0, U1)
        np.testing.assert_almost_equal(V0, V1)
        np.testing.assert_almost_equal(U0, U3,3)
        np.testing.assert_almost_equal(V0, V3,3)

        if plot:
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots(1, 3, sharey=False, figsize=(15.4,3.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.50)
            ax =  flowfield2D_plot(X, Y, U0, V0, bounded=False, xs=xs, ys = ys, ax=axes[0], maxVal=0.5, minVal=0)
            ax =  flowfield2D_plot(X, Y, U1, V1, bounded=False, xs=xs, ys = ys, ax=axes[1], maxVal=0.5, minVal=0)
            ax =  flowfield2D_plot(X, Y, U3, V3, bounded=False, xs=xs, ys = ys, ax=axes[2], maxVal=0.5, minVal=0)
            plt.show()

if __name__ == "__main__":
#     Test().test_CSP_flowrate()
    unittest.main()

