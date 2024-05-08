import numpy as np


def circulation2D(xc, yc, uc, vc, verbose=True):
    """
    Compute the circulation about an arbitrary contour using the trapezoidal rule for integration.
    
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
    INPUTS:
     -  xi (ndarray): Values where velocity is to be recomputed
     -  yi (ndarray): Values where velocity is to be recomputed
     -  u (ndarray): 2D array of x-component velocities.
     -  v (ndarray): 2D array of y-component velocities.
     -  x (ndarray): 2D array of x coordinates.
     -  y (ndarray): 2D array of y coordinates.
     -  method : linear, cubic, nearest
        
    Returns:
        omz (ndarray): 2D array of vorticity.
    """
    import scipy.interpolate as scint
    from matplotlib.tri import Triangulation, LinearTriInterpolator
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

