""" 
Implement a distributed point panel method.

Reference: 
 [1] Anderson Fundamentals-of-aerodynamics-6-Edition, p.292
"""
import numpy as np
from welib.essentials import *
from welib.vortilib.panelcodes.panel_tools import airfoil_params, plot_airfoil
from welib.vortilib.elements.SourcePanel2D import *
from welib.CFD.flows2D import *

def CCSP_panel_solve(XP, YP, Ux, Uy, lift=True, verbose=False):
    r""" 
    Solve the flow about a contiguous (potentially closed-loop) surface using the source panel method.

    INPUTS:
     -XP, YP : coordinates of source points, assumed to go from TE to LE clockwise [m]
     -U_xy: components of the freestream velocity in the x and y direction [m/s]
    OUTPUTS:
     - out: storage for multiple variables, like Cp
    """
    out = {}

    # --- Geometry
    ns = -1 *np.sign(np.sum(XP[:-1]*YP[1:] - XP[1:]*YP[:-1]))
    if ns==0:
        raise Exception('Not a closed countour')
    if verbose :
        print('[INFO] Contour is {}'.format({-1:'counterclockwise', 1:'clockwise'}[ns])) 
    P     = np.column_stack((XP, YP))
    mids  = (P[:-1,:] + P[1:,:]) / 2
    dP    = P[1:,:] - P[:-1,:]
    ds    = np.linalg.norm(dP, axis= 1)
    t_hat = dP / ds[:, np.newaxis]
    n_hat = ns * np.column_stack((-t_hat[:, 1], t_hat[:, 0]))

    # Compute geometric integral
    Vinf =np.sqrt(Ux**2 + Uy**2)
    # --- Right hand side
    rhs = - (Ux * n_hat[:, 0] + Uy * n_hat[:, 1])
    # --- Build matrix
    nS  = len(XP) - 1  # Number of source panels
    SP = P # Source panel points 
    CP = mids
    nCP = len(CP)
    M = np.zeros((nCP,nS))
    for i in range(nCP): # 
        for j in range(nS):
#             if (i == j): 
#                 M[i,j] = 0.5  #  Principal value
#             else:
            u, v = csp_u11(CP[i,:], SP[j,:], SP[j+1,:])
            M[i, j] = (u * n_hat[i, 0] + v * n_hat[i, 1] )
    # --- SOLVE
    sigmas = np.linalg.solve(M, rhs)
    assert(abs(sum(sigmas*ds))<1e-8) 

    # --- Outputs
    # Geometry
    out['x']        = XP
    out['y']        = YP
    out['theta']    = np.arctan2(YP, XP)
    out['ds']       = ds
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['SP']       = P
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





if __name__ == '__main__':
    import matplotlib.pyplot as plt
    np.set_printoptions(linewidth=300, precision=3)
    # --- Geometry
    m  = 151
    m  = 31
    U0 = 1
    R  = 2
    theta   =-np.linspace(0,2*np.pi, m+1)
    xCP = R*np.cos(theta)
    yCP = R*np.sin(theta)
    
    # Panel method
    sigmas, out = CCSP_panel_solve(xCP, yCP, U0, 0)


    # --- Plots
    theta_mid = out['theta_CP']
    Uth_theory = -2.0*U0*np.sin(theta_mid)
    Cp_theory  = 1-(Uth_theory)**2/U0**2
    
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.40)
    axes[0].plot(out['theta_CP'], out['Ut'], label='Ut')
    axes[0].plot(out['theta_CP'], out['Un'], label='Un')
    axes[0].plot(out['theta_CP'], Uth_theory, 'k.', label='Theory')

    axes[0].legend()
    axes[1].plot(out['theta_CP'], out['Cp'], label='Cp')
    axes[1].plot(out['theta_CP'], Cp_theory, 'k.', label='Theory')
    axes[1].legend()
    
    # --- Flow field
    vel = lambda X, Y : ccsp_u(X, Y, out['SP'], out['sigmas'])
    X, Y, U, V =  flowfield2D(vel, xmax=3.5, nx=50, U0x=U0, L=R, rel=True)
    ax =  flowfield2D_plot(X, Y, U, V, ax=axes[2], minVal=0, maxVal=2, bounded=False, rel=True)
    ax.plot(xCP/R, yCP/R, 'k-',lw=3)
    ax.set_title('Source panel, cylinder')
    np.testing.assert_almost_equal(Cp_theory, out['Cp'])

    plt.show()
