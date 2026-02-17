import numpy as np
from welib.vortilib.panelcodes.panel_tools import line_params2
from welib.vortilib.elements.SourcePanel2D import dcsp_u, csp_u11

# NOTE: if velocity is needed, use dcsp_u

def dCSP_solve(SP1, SP2, Uxy=None, fU=None, verbose=False):
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
    Vwall = np.asarray(dcsp_u(CP[:, 0], CP[:, 1], SP1, SP2, sigmas)).T
    Vwall[:,0] += Ux
    Vwall[:,1] += Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # --- Output: Cp
    Vinf2 = Ux**2 + Uy**2 # NOTE: only fine for constant velocity
    out['Cp'] = 1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/Vinf2

    return out
