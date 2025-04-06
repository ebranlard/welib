""" 
Implement a distributed point panel method.

Reference: 
 [1] Anderson Fundamentals-of-aerodynamics-6-Edition, p.292
"""
import numpy as np
from welib.essentials import *
from welib.vortilib.panelcodes.panel_tools import line_params, plot_line, line_params2
from welib.vortilib.elements.SourcePanel2D import *


def CCSP_panel_solve(XP, YP, Uxy=None, fU=None, closed=True, verbose=False):
    r""" 
    Solve the flow about a contiguous (potentially closed-loop) surface using the source panel method.

    INPUTS:
     -XP, YP : coordinates of source points, assumed to go from TE to LE clockwise [m]
     -U_xy: freestream velocity (2 values)
     -fU: freestream velocity function, with interface U,V = fU(X,Y)  (vectorial)
    OUTPUTS:
     - out: storage for multiple variables, like Cp
    """
    out = {}
    # --- Velocity function
    if Uxy is not None:
        fU =lambda X, Y : (X*0+Uxy[0], X*0+Uxy[1])

    # --- Geometry
    ns = -1 *np.sign(np.sum(XP[:-1]*YP[1:] - XP[1:]*YP[:-1]))
    closed_detected=True
    if ns==0:
        if closed:
            print('[WARN] CCSP Not a closed contour')
        closed_detected=False
        ns =1
    if verbose :
        print('[INFO] Contour is {}'.format({-1:'counterclockwise', 1:'clockwise'}[ns])) 
    P     = np.column_stack((XP, YP))
    mids  = (P[:-1,:] + P[1:,:]) / 2
    dP    = P[1:,:] - P[:-1,:]
    ds    = np.linalg.norm(dP, axis= 1)
    t_hat = dP / ds[:, np.newaxis]
    n_hat = ns * np.column_stack((-t_hat[:, 1], t_hat[:, 0]))


    # --- Build matrix
    CP = mids # Positions of control points
    SP = P    # Positions of source panels
    M = np.zeros((len(CP),len(SP)-1))
    for i in range(len(CP)):
        for j in range(len(SP)-1):
#             if (i == j): 
#                 M[i,j] = 0.5  #  Principal value
#             else:
            u, v = csp_u11(CP[i,:], SP[j,:], SP[j+1,:])
            #u, v = ccsp_u([CP[i,0]], [CP[i,1]], SP[j:j+2,:], sigmas=[1])
            M[i,j] = (u*n_hat[i,0] + v*n_hat[i,1])
    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    rhs = -(Ux*n_hat[:,0] + Uy*n_hat[:,1])

    # --- SOLVE
    sigmas = np.linalg.solve(M, rhs)
    if closed:
        try:
            assert(abs(sum(sigmas*ds))<1e-8)
        except:
            print('[WARN] Sum of sigma not zero, contour is likely not closed')

    # --- Outputs
    # Geometry
    out['x']        = XP
    out['y']        = YP
    out['theta']    = np.arctan2(YP, XP)
    out['ds']       = ds
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['SP']       = SP
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    out['rhs']      = rhs
    out['sigmas']   = sigmas

    # --- Output: Velocity at wall
    Vwall = np.asarray(ccsp_u(CP[:, 0], CP[:, 1], SP, sigmas)).T
    Vwall[:,0] += Ux
    Vwall[:,1] += Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # --- Output: Cp
    out['Cp'] = 1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/(Ux**2+Uy**2)

    return sigmas, out

def CSP_panel_solve(SP1, SP2, Uxy=None, fU=None, verbose=False):
    r""" 
    Solve the flow about discontinuous panels using the source panel method.

    INPUTS:
     -SP1: first points of panels
     -SP2: second points of panels
     -U_xy: freestream velocity (2 values)
     -fU: freestream velocity function, with interface U,V = fU(X,Y)  (vectorial)
    OUTPUTS:
     - out: storage for multiple variables, like Cp
    """
    out = {}
    # --- Velocity function
    if Uxy is not None:
        fU =lambda X, Y : (X*0+Uxy[0], X*0+Uxy[1])

    # --- Geometry
    n_hat, t_hat, mids, ds, _, ax =  line_params2(SP1, SP2, plot=False, ntScale=0.3)



    # --- Build matrix
    CP = mids # Positions of control points
    M = np.zeros((len(CP),len(SP1)))
    for i in range(len(CP)):
        for j in range(len(SP1)):
            u, v = csp_u11(CP[i,:], SP1[j], SP2[j,:])
            M[i,j] = (u*n_hat[i,0] + v*n_hat[i,1])
    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    rhs = -(Ux*n_hat[:,0] + Uy*n_hat[:,1])
    #printMat(M)
    #printMat(rhs)

    # --- SOLVE
    sigmas = np.linalg.solve(M, rhs)

    # --- Outputs
    # Geometry
    out['ds']       = ds
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['SP1']      = SP1
    out['SP2']      = SP2
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    out['rhs']      = rhs
    out['sigmas']   = sigmas

    # --- Output: Velocity at wall
    Vwall = np.asarray(csp_u(CP[:, 0], CP[:, 1], SP1, SP2, sigmas)).T
    Vwall[:,0] += Ux
    Vwall[:,1] += Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # --- Output: Cp
    out['Cp'] = 1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/(Ux**2+Uy**2)

    return sigmas, out




if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from welib.CFD.flows2D import *    
    np.set_printoptions(linewidth=300, precision=3)
    method='CCSP'
    method='CSP'
    # --- Geometry
    m  = 151
    m  = 51
    U0 = 1
    R  = 2
    theta   =-np.linspace(0,2*np.pi, m+1)
    XP = R*np.cos(theta)
    YP = R*np.sin(theta)

    
    # --- Panel method
    if method=='CCSP':
        sigmas, out = CCSP_panel_solve(XP, YP, Uxy=(U0, 0))
    else:
        SP1 = np.column_stack([XP[0:-1], YP[0:-1]])
        SP2 = np.column_stack([XP[1:]  , YP[1:]])
        sigmas, out = CSP_panel_solve(SP1, SP2, Uxy=(U0, 0))

    # --- Theory
    theta_mid = out['theta_CP']
    Ut_theory = -2.0*U0*np.sin(theta_mid)
    Cp_theory  = 1-(Ut_theory)**2/U0**2
    
    # --- Plots   
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.40)
    ax = axes[0]

    ax.plot(out['theta_CP'], out['Un'], label='Un')
    ax.plot(out['theta_CP'],-out['Ut'], label='Ut')
    ax.plot(out['theta_CP'], Ut_theory, 'k.', label='Theory')

    ax.legend()
    ax = axes[1]
    ax.plot(out['theta_CP'], out['Cp'], label='Cp')
    ax.plot(out['theta_CP'], Cp_theory, 'k.', label='Theory')
    ax.legend()
    
    # --- Flow field
    if method=='CCSP':
        vel = lambda X, Y : ccsp_u(X, Y, out['SP'], out['sigmas'])
    else:
        vel = lambda X, Y : csp_u(X, Y, SP1, SP2, out['sigmas'])
    X, Y, U, V =  flowfield2D(vel, xmax=3.5, nx=50, U0x=U0, L=R, rel=True)
    ax =  flowfield2D_plot(X, Y, U, V, ax=axes[2], minVal=0, maxVal=2, bounded=False, rel=True)
    ax.plot(XP/R, YP/R, 'k-',lw=3)
    ax.set_title('Source panel, cylinder '+method)
    np.testing.assert_almost_equal(Cp_theory, out['Cp'])

    plt.show()
