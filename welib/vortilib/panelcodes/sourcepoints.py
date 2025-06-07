""" 
Implement a source point panel method.

"""
import numpy as np
from welib.essentials import *
from welib.vortilib.elements.SourcePoint import *
from welib.vortilib.panelcodes.panel_tools import line_params, plot_line, line_params2


# --------------------------------------------------------------------------------}
# --- 2d point source 
# --------------------------------------------------------------------------------{
def source_flow(DX, DY, Sigma=1): 
    """ 
    Induced velocity by one 2D source point on multiple Control Points (CP)
    INPUTS:
    - CP: position of control point
    - Ps: position of source

    EXAMPLE:
        # [Ui]= sp2d_u(1, 0, 2*np.pi)
        # [Ui]= sp2d_u(0, 1, 2*np.pi)
    """
    r2 = DX ** 2 + DY ** 2
    rX = DX
    rY = DY
    bOK=r2>1e-5
    tU = np.zeros_like(DX)
    V = np.zeros_like(DX)
    U = np.zeros_like(DX)
    V = np.zeros_like(DX)
    U[bOK] = Sigma/(2*np.pi) * rX[bOK] / r2[bOK]   # [1] Eq. 32.3
    V[bOK] = Sigma/(2*np.pi) * rY[bOK] / r2[bOK] 
#     U[~bOK] = 0.5*Sigma* rX[~bOK] / r2[~bOK]
#     V[~bOK] = 0.5*Sigma* rY[~bOK] / r2[~bOK]
    return U, V

def CSP_velocity(X, Y, SP, Sigmas, debug=False):
    nS = len(Sigmas)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(nS):
        sp = SP[i]
        us, vs = source_flow(X-sp[0], Y-sp[1], Sigma=Sigmas[i])
        U+=us
        V+=vs
    return U, V
sp_u = CSP_velocity


def CSP_panel_solve(XP, YP, Uxy=None, fU=None, verbose=False, offset=0.1):
    r""" 
    Solve the flow about a contiguous (potentially closed-loop) surface using the point source panel method.

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
    n_hat, t_hat, mids, ds, _, ax =  line_params(XP, YP, plot=False, ntScale=0.3)


    # --- Build matrix
    nCP = len(XP)-1 # Number of control points
    nS  = len(XP)-1 # Number of source panels
    CP = mids # Positions of control points
    SP = mids # Positions of sources
    if offset >0:
        SP = mids - offset*(n_hat.T*ds).T #-0.05*t_hat
    M = np.zeros((nCP,nS))
    for i in range(nCP):
        for j in range(nS):
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

    # --- SOLVE
    Sigmas = np.linalg.solve(M, rhs)
    sigmas = Sigmas/ds # strength per unit length
    assert(abs(sum(Sigmas))<1e-8)

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
    out['Sigmas']   = Sigmas

    # --- Output: Velocity at wall
    Vwall = np.zeros_like(mids)
    Vwall[:,0], Vwall[:,1] = CSP_velocity(CP[:,0], CP[:,1], SP, Sigmas)
    Vwall[:,0] += Ux
    Vwall[:,1] += Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # --- Output: Cp
    out['Cp'] = 1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/(Ux**2+Uy**2)
    return Sigmas, out


def SP_panel_solve(SP1, SP2, Uxy=None, fU=None, offset=0.1, verbose=False):
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



    # --- Build matrix
    CP = mids # Positions of control points
    SP = mids # Positions of sources
    if offset >0:
        SP = mids - offset*(n_hat.T*ds).T #-0.05*t_hat
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

    # --- SOLVE
    Sigmas = np.linalg.solve(M, rhs)
    sigmas = Sigmas/ds # strength per unit length
    assert(abs(sum(Sigmas))<1e-8)

    # --- Outputs
    # Geometry
    out['x']        = SP[:,0]
    out['y']        = SP[:,1]
    out['theta']    = np.arctan2(SP[:,0], SP[:,1])
    out['ds']       = ds
    out['n']        = n_hat
    out['t']        = t_hat
    out['CP']       = CP
    out['SP']       = SP
    out['SP1']      = SP1
    out['SP2']      = SP2
    out['theta_CP'] = np.arctan2(CP[:,1], CP[:,0])
    out['rhs']      = rhs
    out['sigmas']   = sigmas
    out['Sigmas']   = Sigmas

    # --- Output: Velocity at wall
    Vwall = np.zeros_like(mids)
    Vwall[:,0], Vwall[:,1] = CSP_velocity(CP[:,0], CP[:,1], SP, Sigmas)
    Vwall[:,0] += Ux
    Vwall[:,1] += Uy
    out['Vwall'] = Vwall
    out['Un'] = Vwall[:,0]*n_hat[:,0] + Vwall[:,1]*n_hat[:,1]
    out['Ut'] = Vwall[:,0]*t_hat[:,0] + Vwall[:,1]*t_hat[:,1]
    # --- Output: Cp
    out['Cp'] = 1 - (Vwall[:,0]**2 + Vwall[:,1]**2)/(Ux**2+Uy**2)
    return Sigmas, out




if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from welib.CFD.flows2D import *    
    method='CSP'
    method='SP'
    case = 'cylinder'
#     case = 'ellipse_lift'
    offset = 2.0
    if case=='cylinder':
        m = 30
        U0    = 1
        R=2
        theta   =-np.linspace(0,2*np.pi, m+1)
        XP = R*np.cos(theta)
        YP = R*np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        Ut =  2.0*U0*np.sin(theta_mid) # gamma_exact = Ut at CP
        maxVal =  2
    elif case=='ellipse_lift':
        m     = 130
        U0    = 1
        R     = 1    # half major
        ratio = 0.5
        theta   =-np.linspace(0,2*np.pi, m+1)
        abyr= np.sqrt((1.-ratio)/(1.+ratio))# ! see Lewis p 50
        XP = R      * np.cos(theta)
        YP = R*ratio* np.sin(theta)
        dtheta=theta[1]-theta[0]
        theta_mid   = theta[:-1]+dtheta/2
        Ut =  (2.0*U0*np.sin(theta_mid))/np.sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*np.cos(2.*theta_mid) )  # See Lewis p50
        maxVal =  1.5

    # --- Panel method
    if method=='CSP':
       Sigmas, out = CSP_panel_solve(XP, YP, (U0, 0), verbose=True, offset=offset)
    else:
        SP1 = np.column_stack([XP[0:-1], YP[0:-1]])
        SP2 = np.column_stack([XP[1:]  , YP[1:]])
        sigmas, out = SP_panel_solve(SP1, SP2, Uxy=(U0, 0), offset=offset)
   
    # --- Theory
    # Note for ellipst, we can't use "theta" from the CP
    Ut_theory = Ut
    Cp_theory = 1-Ut**2/U0**2
    
    # --- Plots   
    fig, axes = plt.subplots(1, 4, figsize=(10, 3.5))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.40)

    ax = plot_line(XP, YP, Uwall=out['Vwall'], ntScale=0.5, UScale=0.8, ax=axes[0])
    ax.plot(out['SP'][:,0], out['SP'][:,1], '+')

    ax=axes[1]
    ax.plot(out['theta_CP'], out['Un'], label='Un')
    ax.plot(out['theta_CP'], out['Ut'], label='Ut')
    ax.plot(out['theta_CP'], Ut_theory, 'k.', label='Theory')

    ax.legend()
    ax=axes[2]
    ax.plot(out['theta_CP'], out['Cp'], label='Cp')
    ax.plot(out['theta_CP'], Cp_theory, 'k.', label='Theory')
    ax.legend()
    
    # --- Flow field
    vel = lambda X, Y : CSP_velocity(X, Y, out['SP'], out['Sigmas'])
    X, Y, U, V =  flowfield2D(vel, xmax=3.5, nx=50, U0x=U0, L=R, rel=True)
    if case=='cylinder':
        RR =np.sqrt(X**2+Y**2)
        U[RR<1]=0
        V[RR<1]=0
    ax =  flowfield2D_plot(X, Y, U, V, ax=axes[3], minVal=0, maxVal=maxVal, bounded=False, rel=True)
    ax.plot(XP/R, YP/R, 'k-',lw=3)
    #ax.set_title('Source panel, cylinder')

    print('MaxError ',np.max(np.abs(Cp_theory-out['Cp'])))
#     print('')
#     np.testing.assert_almost_equal(Cp_theory, out['Cp'], 2)



    plt.show()
