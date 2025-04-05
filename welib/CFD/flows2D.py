import numpy as np


def circulation2D(xc, yc, uc, vc, verbose=True):
    r"""
    Compute the circulation about an arbitrary contour using the trapezoidal rule for integration.

    Circulation = \int  u · that ds
    
    INPUTS:
    - xc: 1D array containing the x-coordinates of the contour nodes.
    - yc: 1D array containing the y-coordinates of the contour nodes.
    - uc: 1D array containing the x-components of the velocity at each node.
    - vc: 1D array containing the y-components of the velocity at each node.
    
    OUTPUTS:
    - Circulation around the contour.
    """
    # Check if the last point is the same as the first point
    if np.abs(xc[0]-xc[-1])>1e-9 or np.abs(yc[0]-yc[-1])>1e-9:
        if verbose: 
            print('[WARN] circulation2D: Contour was not complete > appending first point.')
        # Add the first point as the last point
        xc = np.append(xc, xc[0])
        yc = np.append(yc, yc[0])
        uc = np.append(uc, uc[0])
        vc = np.append(vc, vc[0])

    # Compute the differential element dl
    dlx = np.diff(xc)
    dly = np.diff(yc)
    # Compute the average velocity at each node
    u_mid = 0.5 * (uc[:-1] + uc[1:])
    v_mid = 0.5 * (vc[:-1] + vc[1:])
    
    # Compute the scalar product of velocity and dl
    scalar_product = u_mid * dlx + v_mid * dly
    circulation = np.sum(scalar_product)
    return circulation

def flowrate2D(xc, yc, uc, vc, verbose=False, ns=1):
    r"""
    Compute the net flow rate across a 2D line

    Flow Rate = \int  u · nhat ds
                \iint  nabla\cdot u  dA    < for closed contours

    Note:
        n ds = [-dy/ds , dx/ds] ds = [-dy, dx]

    INPUTS:
     - xc, yc : coordinates of line points (may form a loop if last==first)
     - uc, vc : velocity components at each point
     - verbose: if True, prints intermediate info
    OUTPUTS:
     - Q: Net flow rate 
    """
    # --- Sign
    #ns = -1 *np.sign(np.sum(xc[:-1]*yc[1:] - xc[1:]*yc[:-1]))
    #if verbose :
    #    print('[INFO] Contour is {}'.format({-1:'counterclockwise', 1:'clockwise'}[ns])) 
    # --- Geometry
    dlx = np.diff(xc)
    dly = np.diff(yc)
    # Compute the average velocity at each node
    u_mid = 0.5 * (uc[:-1] + uc[1:])
    v_mid = 0.5 * (vc[:-1] + vc[1:])
    # Compute the scalar product of velocity and (n ds) = [-dy/ds , dx/ds] ds = [-dy, dx]
    scalar_product = ns * (- u_mid * dly + v_mid * dlx)
    Q = np.sum(scalar_product)
    return Q


def vorticity2D(u, v, x, y):
    """
    Compute the vorticity (ωz) based on u and v velocities.
    
    Parameters:
        u (ndarray): 2D array of x-component velocities.
        v (ndarray): 2D array of y-component velocities.
        x (ndarray): 2D array of x coordinates.
        y (ndarray): 2D array of y coordinates.
        
    Returns:
        omz (ndarray): 2D array of vorticity.
    """

    x_flat = np.sort(np.unique(x.flatten()))
    y_flat = np.sort(np.unique(y.flatten()))

    # Compute the derivatives of u and v with respect to y and x, respectively
    # Option 1
    du_dy, du_dx = np.gradient(u, y_flat, x_flat, axis=(0,1))
    dv_dy, dv_dx = np.gradient(v, y_flat, x_flat, axis=(0,1))
    # Compute the spatial derivatives of u with respect to y and v with respect to x
    # Option 2
    #du_dy = np.gradient(u, y_flat, axis=0)
    #dv_dx = np.gradient(v, x_flat, axis=1)
    
    # Compute vorticity (ωz) using the formula: ωz = dv/dx - du/dy 
    omz = dv_dx - du_dy 
    
    return omz

def flow_interp2D(xi, yi, u, v, x, y, method='linear', algo='griddata'):
    """
    The flow may be known on a grid or not.
    INPUTS:
     -  xi (ndarray): Values where velocity is to be recomputed
     -  yi (ndarray): Values where velocity is to be recomputed
     -  u (ndarray): 1D or 2D array of x-component velocities.
     -  v (ndarray): 1D or 2D array of y-component velocities.
     -  x (ndarray): 1D or 2D array of x coordinates.
     -  y (ndarray): 1D or 2D array of y coordinates.
     -  method : linear, cubic, nearest
        
    Returns:
        omz (ndarray): 2D array of vorticity.
    """
    import scipy.interpolate as scint
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    x = np.asarray(x)
    y = np.asarray(y)
    u = np.asarray(u)
    v = np.asarray(v)
    if algo=='griddata':
        x=x.flatten()
        y=y.flatten()
        u=u.flatten()
        v=v.flatten()
        Ui = scint.griddata((x, y), u, (xi,yi), method=method)
        Vi = scint.griddata((x, y), v, (xi,yi), method=method)
    elif algo=='TriInterpolator':
        x=x.flatten()
        y=y.flatten()
        u=u.flatten()
        v=v.flatten()
        tri = Triangulation(x, y)
        if method=='linear':
            fu = LinearTriInterpolator(tri, u)
            fv = LinearTriInterpolator(tri, v)
        elif method=='cubic':
            fu = CubicTriInterpolator(tri, u)
            fv = CubicTriInterpolator(tri, v)
        else:
            raise NotImplementedError()
        Ui = fu(xi, yi)
        Vi = fv(xi, yi)
    else:
        raise NotImplementedError()
    return Ui, Vi


def flowfield2D(function, xmax=1, ymax=None, xmin=None, ymin=None, nx = 50, ny = None, U0x=0, U0y=0, Vref=None, L=1, rel=False):
    """ Evaluate a function to get a velocity field on a grid
    INPUTS:
      - function: function with interface U,V = function(X,Y)
    """
    # --- Default is a squared grid
    if ymax is None:
        ymax = xmax
    if xmin is None:
        xmin = -xmax
    if ymin is None:
        ymin = -ymax
    if ny is None:
        ny = nx+1
    vx = np.linspace(xmin, xmax, nx)
    vy = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(vx, vy)
    U, V = function(X, Y)
    U += U0x
    V += U0y
    if rel:
        # Relative 
        if Vref is None:
            Vref = np.sqrt(U0x+U0y)
        U = U/Vref
        V = V/Vref
        X= X/L
        Y= Y/L
    return X, Y, U, V

def flowfield2D_plot(
        X, Y, U, V,
        ax = None,
        minVal = None, maxVal = None, nLevels=11, bounded=False,
        ctOpts=None,
        stOpts=None, xs=None, ys=None,
        speed = True,
        xlabel =None,
        ylabel =None,
        clabel =None,
        rel = False,
        ):
    """ """
    import matplotlib.pyplot as plt
    # --- Default options
    if ctOpts is None:
        ctOpts={}
    if stOpts is None:
        stOpts={}
    for k, v in {'linewidth':0.7, 'density':30, 'arrowstyle':'->'}.items():
        if k not in stOpts:
            stOpts.update({k:v})
    if clabel is None:
        clabel = r'$\sqrt{u^2+v^2}/U_\text{ref}$' if rel else r'$\sqrt{u^2+v^2}$'
    if xlabel is None:
        xlabel = r'$x/L$ [-]' if rel else r'$x [m]$'
    if ylabel is None:
        ylabel = r'$y/L$ [-]' if rel else r'$y [m]$'

    # --- Axes
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    else:
        fig = ax.figure

    # --- Ranges
    xmin = np.min(X.flatten())
    xmax = np.max(X.flatten())
    ymin = np.min(Y.flatten())
    ymax = np.max(Y.flatten())
    vx = X[0,:]
    vy = Y[:,0]
    # --- 
    Speed = np.sqrt((U**2+V**2))
    if minVal is None:
        minVal = np.min(Speed.flatten())
    if maxVal is None:
        maxVal = np.max(Speed.flatten())
    if bounded:
        Speed[Speed<minVal] = minVal
        Speed[Speed>maxVal] = maxVal

    #  --- Streamlines
    if xs is None and ys is None:
        ys = np.linspace(ymin, ymax, 15)
        xs = ys*0 -xmax
    elif xs is None:
        xs = ys*0 -xmax
    elif ys is None:
        if not hasattr(xs,'__len__'):
            xs = [xs]*15
        ys = np.linspace(ymin, ymax, len(xs))
    start = np.array([xs, ys])
    start = np.array([[xi, yi ] for xi,yi in start.T if xmin<=xi<=xmax and ymin<=yi<=ymax]).T
    if len(start)==0:
        raise Exception('Streamlines are not within domain')

    # --- Plot
    im = ax.contourf(X, Y, Speed, levels=np.linspace(minVal, maxVal, nLevels), vmin=minVal, vmax=maxVal, **ctOpts)
    cb = fig.colorbar(im)
    cb.set_label(clabel)
    sp = ax.streamplot(vx, vy, U, V, color='k', start_points=start.T, **stOpts)
    #sp = ax.streamplot(vx, vy, U, V, color='k', linewidth=0.7,density=30,arrowstyle='-')
    #qv = streamQuiver(ax, sp, spacing=1, offset=1.0, scale=40, angles='xy')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect('equal','box')
    ax.tick_params(direction='in', top=True, right=True)
    return ax



if __name__ == '__main__':
    theta = np.linspace(0, 2*np.pi, 250, endpoint=True)
    xc = np.cos(theta)
    yc = np.sin(theta)
    uc = xc  # radial field: u = x
    vc = yc  #               v = y
    Q = flowrate2D(xc, yc, uc, vc, verbose=True)
    print(f"Flow rate: {Q:.4f} {2*np.pi:.4f}")
