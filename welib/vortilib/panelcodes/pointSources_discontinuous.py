""" 
Implement a source point panel method.

"""
import numpy as np
from welib.essentials import *
#from welib.vortilib.elements.SourcePoint import *
from welib.vortilib.panelcodes.panel_tools import line_params, plot_line, line_params2
from welib.vortilib.panelcodes.pointSources import source_flow, PS_velocity

def dPS_velocity(X, Y, SP1, SP2, Sigmas, debug=False):
    mid = (SP1 + SP2)/2
    return PS_velocity(X, Y, mid, Sigmas)


def dPS_solve(SP1, SP2, Uxy=None, fU=None, offset=0.1, verbose=False):
    r""" 
    Solve the flow about discontinuous panels using the point source panel method.

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
    
    # --- Panels / Control points  
    CP = mids # Positions of control points
    SP = mids # Positions of sources
    if offset >0:
        SP = mids - offset*(n_hat.T*ds).T #-0.05*t_hat

    # --- Build system matrix
    M = np.zeros((len(CP),len(SP)))
    for i in range(len(CP)):
        for j in range(len(SP)):
            if offset==0 and i==j:
                M[i,i] = 0.5 *ds[i]
            else:
               rj = SP[j,:] # 
               ri = CP[i,:] # 
               dr = ri-rj
               u, v = source_flow(dr[0], dr[1], Sigma=1)
               M[i,j] = (u*n_hat[i,0] + v*n_hat[i,1])

    # --- Right hand side
    Ux, Uy = fU(CP[:,0], CP[:,1])
    rhs = -(Ux*n_hat[:,0] + Uy*n_hat[:,1])

    # --- Solve
    Sigmas = np.linalg.solve(M, rhs)
    sigmas = Sigmas/ds # strength per unit length
    #assert(abs(sum(Sigmas))<1e-8)
    
    # Tangential velocities on panel
    Vwall = np.zeros_like(mids)
    Vwall[:,0], Vwall[:,1] = PS_velocity(CP[:,0], CP[:,1], SP, Sigmas)
    Vwall[:,0] += Ux
    Vwall[:,1] += Uy

    # Pressure coefficient
    Vinf2 = Ux**2 + Uy**2 # NOTE: only fine for constant velocity        
    Cp = 1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/Vinf2

    # --- Outputs
    # Geometry
    out['x']        = SP[:,0]
    out['y']        = SP[:,1]
    out['theta']    = np.arctan2(SP[:,0], SP[:,1])
    out['ds']       = ds        # Panel lengths
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['SP']       = SP
    out['SP1']      = SP1
    out['SP2']      = SP2
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    # Output: System
    out['rhs']      = rhs       # Right hand side
    out['M']        = M         # System matrix
    # Output: Solution
    out['sigmas']   = sigmas
    out['Sigmas']   = Sigmas

    # --- Output: Velocity at wall

    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # Output: Cp
    out['Cp'] = Cp
    return out


if __name__ == '__main__':
    # See continuous example.
    plt.show()
