""" 
Set of tools useful for 2D panel methods

"""

import numpy as np


# --------------------------------------------------------------------------------
# --- Geometry 
# --------------------------------------------------------------------------------
def compute_curvature(X, Y, method='Menger'):
    """ """
    nP = len(X)
    curv     = np.zeros(nP-1)
    tangents = np.zeros((nP-1,2))
    ds       = np.zeros(nP-1)
    for i in range(nP-1):
        P1 = np.array((X[i],Y[i]))
        P2 = np.array((X[i+1],Y[i+1]))
        dP = P2-P1
        ds[i]       = np.linalg.norm(dP)
        tangents[i] = dP/ds[i]
    if method=='Menger':
        # --- Menger
        # Compute Menger curvature (circle passing through three points)
        X = np.concatenate(([X[-1]],X,[X[0]]))
        Y = np.concatenate(([Y[-1]],Y,[Y[0]]))
        for i in range(1,nP):
            P1 = np.array((X[i-1],Y[i-1]))
            P2 = np.array((X[i]  ,Y[i]))
            P3 = np.array((X[i+1],Y[i+1]))
            L1= np.linalg.norm(P2-P1)
            L2= np.linalg.norm(P3-P2)
            L3= np.linalg.norm(P1-P3)
            area = 0.5*((P2[0] - P1[0]) * (P3[1] -P1[1]) - (P2[1] - P1[1]) * (P3[0] -P1[0]))
            curv[i-1] = 4*area/(L1*L2*L3)
    elif method=='Lewis':
        # --- Lewis
        # WATCH OUT THIS MOSTLY WORK WITH LEWIS "PHI" angle going from LE to TE by suction
        # 1/rm = DeltaBeta_m/ds [Lewis 1.31, p 25]
        # cds[:] = -2.51327 # m=5
        # cds[:] = -0.66138793 # m=19
        m=len(X)-1
        slope = np.arctan2(tangents[:,1], tangents[:,0])
        b = slope>np.pi/2
        slope[b] = slope[b]-2*np.pi
        cds = slope*0
        cds[0]     = slope[1]-slope[-1] -2*np.pi
        cds[1:m-1] = slope[2:]-slope[0:-2]
        cds[m-1]   = slope[0]-slope[-2] -2*np.pi
        curv = cds / ds /2
    elif method=='zero':
        pass
    else:
        raise NotImplementedError('Method:`{}`'.format(method))
    bNaN = np.isnan(curv)
    if sum(bNaN)>0:
        print('[WARN] panel_tools: compute_curvature, {} nan encountered'.format(sum(bNaN)))
        curv[bNaN]=0
    return curv

# --------------------------------------------------------------------------------}
# --- Line / airfoil
# --------------------------------------------------------------------------------{
def plot_airfoil(*args, **kwargs):
    return plot_line(*args, **kwargs)

def airfoil_params(*args, **kwargs):
    return line_params(*args, **kwargs)

def line_params(X, Y, plot=False, ntScale=0.3, curv_method='Menger', verbose=False):
    """ 
    Compute normals, tangents, midpoint and ds for an airfoil
    The coordinates are assumed to go from lower TE to upper TE clockwise
    INPUTS: 
      - X: array of x coordinates, size n
      - Y: array of y coordinates, size n
    OUTPUTS:
      - normals: array of normal vectors, size nx2
      - tangents: array of normal vectors, size nx2
      - mids   : array of mid point coordinats, size nx2
      - ds   : array of panel length, size n
      - ax   : axis if a plot is generated
    """

    # --- Detect clockwise if close countour only
    ns = -1 *np.sign(np.sum(X[:-1]*Y[1:] - X[1:]*Y[:-1]))
    if ns==0:
        print('[WARN] CCSP Not a closed contour')
        ns =1
    if verbose :
        print('[INFO] Contour is {}'.format({-1:'counterclockwise', 1:'clockwise'}[ns])) 

    # --- Geometry
    P     = np.column_stack((X, Y))
    mids  = (P[:-1,:] + P[1:,:]) / 2
    dP    = P[1:,:] - P[:-1,:]
    ds    = np.linalg.norm(dP, axis= 1)
    t_hat = dP / ds[:, np.newaxis]
    n_hat = ns * np.column_stack((-t_hat[:, 1], t_hat[:, 0]))
    if any(ds<1e-8): 
        raise Exception('Some Panels very small')

    curv = compute_curvature(X, Y, method=curv_method)


    if plot:
        maxDs = np.max(ds)
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(X, Y)
        scale=maxDs*ntScale
        for Pmid,t,n in zip(mids, t_hat, n_hat):
            ax.plot(  Pmid[0]+np.array([0, n[0]])*scale, Pmid[1]+np.array([0, n[1]])*scale, 'k')
            ax.plot(  Pmid[0]+np.array([0, t[0]])*scale, Pmid[1]+np.array([0, t[1]])*scale, 'k')
        ax.plot(X[0], Y[0], 's')
        ax.plot(X[1], Y[1], 'o')
        ax.set_aspect('equal', 'box')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
    else:
        ax = None
    return n_hat, t_hat, mids, ds, curv, ax


    

def plot_line(X, Y, Uwall=None, ntScale=0.1, UScale=0.1, nt=True, ax=None):
    import matplotlib.pyplot as plt

    normals, tangents, mids, ds, _, _ = line_params(X, Y, plot=False, curv_method='zero')

    if ax is None:
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    else:
        fig = ax.figure

    maxDs = np.max(ds)

    # --- Plot
    ax.plot(X, Y)

    if Uwall is not None:
        scale=maxDs*UScale
        for i,(Pmid,t,n) in enumerate(zip(mids, tangents, normals)):
            ax.plot(  Pmid[0]+np.array([0, Uwall[i,0]])*scale, Pmid[1]+np.array([0, Uwall[i,1]])*scale, 'r')


    if nt is not None:
        scale=maxDs*ntScale
        for Pmid,t,n in zip(mids, tangents, normals):
            ax.plot(  Pmid[0]+np.array([0, n[0]])*scale, Pmid[1]+np.array([0, n[1]])*scale, 'k')
            ax.plot(  Pmid[0]+np.array([0, t[0]])*scale, Pmid[1]+np.array([0, t[1]])*scale, 'k')
    ax.plot(X[0], Y[0], 's')
    ax.plot(X[1], Y[1], 'o')
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')

    return ax

def line_params2(P1, P2, plot=False, ntScale=0.3, curv_method='Menger', verbose=False):
    mids  = (P1 + P2) / 2
    dP    = P2- P1
    ds    = np.linalg.norm(dP, axis= 1)
    t_hat = dP / ds[:, np.newaxis]
    n_hat = np.column_stack((-t_hat[:, 1], t_hat[:, 0]))
    if any(ds<1e-8): 
        raise Exception('Some Panels very small')
    return n_hat, t_hat, mids, ds, None, None

def plot_line2(P1, P2, Uwall=None, ntScale=0.1, UScale=0.1, nt=True, ax=None):
    from welib.tools.colors import python_colors

    normals, tangents, mids, ds, _, _ = line_params2(P1,P2)

    maxDs = np.max(ds)

    # --- Plot
    for j in range(len(P1)):
        ax.plot([P1[j,0], P2[j,0]], [P1[j,1], P2[j,1]],'-', c= python_colors(1))

    if Uwall is not None:
        scale=maxDs*UScale
        for i,(Pmid,t,n) in enumerate(zip(mids, tangents, normals)):
            ax.plot(  Pmid[0]+np.array([0, Uwall[i,0]])*scale, Pmid[1]+np.array([0, Uwall[i,1]])*scale, 'r')

    if nt is not None:
        scale=maxDs*ntScale
        for Pmid,t,n in zip(mids, tangents, normals):
            ax.plot(  Pmid[0]+np.array([0, n[0]])*scale, Pmid[1]+np.array([0, n[1]])*scale, 'k')
            ax.plot(  Pmid[0]+np.array([0, t[0]])*scale, Pmid[1]+np.array([0, t[1]])*scale, 'k')
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    return ax

