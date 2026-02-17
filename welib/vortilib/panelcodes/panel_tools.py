""" 
Set of tools useful for 2D panel methods

"""

import numpy as np
import matplotlib.pyplot as plt


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
        if np.abs(X[0]-X[-1])<1e-12 or np.abs(Y[0]-Y[-1])<1e-12:
            # Contour is closed
            X = np.concatenate(([X[-2]],X,[X[1]]))
            Y = np.concatenate(([Y[-2]],Y,[Y[1]]))
        else:
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
            if L1*L2*L3 == 0:
                import pdb; pdb.set_trace()
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

def panel_geometry(XP, YP, closed_expected=True, force_clockwise=True):
    ns_in  = -1 * np.sign(np.sum(XP[:-1]*YP[1:] - XP[1:]*YP[:-1]))
    if ns_in==0:
        if closed_expected:
            print('[WARN] Not a closed contour, make sure order makes sense')
        ns = 1
    elif ns_in==-1:
        if force_clockwise:
            print('[INFO] Making contour clockwise')
            XP = XP[::-1]
            YP = YP[::-1]

    ns     = -1 * np.sign(np.sum(XP[:-1]*YP[1:] - XP[1:]*YP[:-1]))
    if ns==0:
        ns=1 # Not a closed countour
    PP     = np.column_stack((XP, YP))
    mids   = (PP[:-1,:] + PP[1:,:]) / 2
    dP     = PP[1:,:] - PP[:-1,:]
    ds     = np.linalg.norm(dP, axis=1)
    t_hat  = dP / ds[:, None]
    n_hat  = ns * np.column_stack((-t_hat[:,1], t_hat[:,0])) # ns Ensures normals are outwards
    phi    = np.arctan2(dP[:,1], dP[:,0])
    phi[phi<0] += 2*np.pi
    return PP, mids, dP, ds, t_hat, n_hat, phi, ns

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


# --------------------------------------------------------------------------------}
# --- Airfoils 
# --------------------------------------------------------------------------------{
def points_inside(XP, YP, XX, YY):
    from matplotlib import path
    AF     = np.vstack((XP.T,YP.T)).T
    afPath = path.Path(AF)
    points = np.column_stack((XX.ravel(), YY.ravel()))
    inside = afPath.contains_points(points).reshape(XX.shape)
    return inside

# --------------------------------------------------------------------------------}
# --- Line / airfoil
# --------------------------------------------------------------------------------{
def plot_airfoil(*args, **kwargs):
    return plot_line(*args, **kwargs)

def airfoil_params(*args, **kwargs):
    return line_params(*args, **kwargs)



    

def plot_line(X, Y, Uwall=None, ntScale=0.1, UScale=0.1, nt=True, ax=None):

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


# --------------------------------------------------------------------------------}
# --- Airfoils Cp
# --------------------------------------------------------------------------------{
def plot_pressure_force_bars(CP, Cp, n_hat, scale=0.1, ax=None):
    """ Plot Pressure force as blue/red normal bars  """
    if ax is None:
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    scale = 0.1
    Cp_abs = np.abs(Cp * scale)
    for i in range(len(CP)):
        x_bar = CP[i,0] + np.array([ 0,  Cp_abs[i]*n_hat[i,0]])
        y_bar = CP[i,1] + np.array([ 0,  Cp_abs[i]*n_hat[i,1]])
        sty='b-' if (Cp[i] > 0) else 'r-'
        ax.plot(x_bar, y_bar, sty)
    #ax.fill(XP,YP,'k')
    ax.set_xlabel('x/c [-]')
    ax.set_ylabel('y/c [-]')
    #ax.legend()
    ax.set_aspect('equal')                                               # Set aspect ratio equal
    return ax


def plot_Cp(x, Cp, ax=None, Cp_ref=None, simple=True, label=None, sty='--', x_ref=None):
    # --- Plot Cp
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    midIndS = int(np.floor(len(Cp)/2))                                          # Airfoil middle index for VPM data

    if Cp_ref is not None:
        if x_ref is None:
            x_ref = x
        ax.plot(x_ref, Cp_ref[:],  'k-', label='Reference Cp')
    if simple:
        ax.plot(x, Cp, sty, label=label)
    else:
        ax.plot(x[midIndS+1:len(x)],Cp[midIndS+1:len(CP)], 'ks', markerfacecolor='b', label=label+' Upper' if label is not None else 'Upper')
        ax.plot(x[0:midIndS       ],        Cp[0:midIndS], 'ks', markerfacecolor='r', label=label+' Lower' if label is not None else 'Lower')


    # ax.set_xlim([0,1])
    ax.set_xlabel('x/c [-]') 
    ax.set_ylabel('Cp [-]')
    ax.legend()
    ax.invert_yaxis()
    return ax
