import numpy as np


def circulation(xc, yc, zc, uc, vc, wc, verbose=True):
    r"""
    Compute the circulation about an arbitrary contour using the trapezoidal rule for integration.

    Circulation = \int  u · that ds
    
    INPUTS:
    - xc: 1D array containing the x-coordinates of the contour nodes.
    - yc: 1D array containing the y-coordinates of the contour nodes.
    - zc: 1D array containing the z-coordinates of the contour nodes.
    - uc: 1D array containing the x-components of the velocity at each node.
    - vc: 1D array containing the y-components of the velocity at each node.
    - wc: 1D array containing the z-components of the velocity at each node.
    
    OUTPUTS:
    - Circulation around the contour.
    """
    # Check if the last point is the same as the first point
    if np.abs(xc[0]-xc[-1])>1e-9 or np.abs(yc[0]-yc[-1])>1e-9 or np.abs(zc[0]-zc[-1])>1e-9: 
        if verbose: 
            print('[WARN] circulation2D: Contour was not complete > appending first point.')
        # Add the first point as the last point
        xc = np.append(xc, xc[0])
        yc = np.append(yc, yc[0])
        zc = np.append(zc, zc[0])
        uc = np.append(uc, uc[0])
        vc = np.append(vc, vc[0])
        wc = np.append(wc, wc[0])

    # Compute the differential element dl
    dlx = np.diff(xc)
    dly = np.diff(yc)
    dlz = np.diff(zc)
    # Compute the average velocity at each node
    u_mid = 0.5 * (uc[:-1] + uc[1:])
    v_mid = 0.5 * (vc[:-1] + vc[1:])
    w_mid = 0.5 * (wc[:-1] + wc[1:])
    
    # Compute the scalar product of velocity and dl
    scalar_product = u_mid * dlx + v_mid * dly + w_mid * dlz
    circulation = np.sum(scalar_product)
    return circulation


def flow_interp(xi, yi, zi, u, v, w, x, y, z, method='linear', algo='griddata'):
    """
    The flow may be known on a grid or not.
    INPUTS:
     -  xi (ndarray): Values where velocity is to be recomputed
     -  yi (ndarray): Values where velocity is to be recomputed
     -  zi (ndarray): Values where velocity is to be recomputed
     -  u (ndarray): 1D or 2D array of x-component velocities.
     -  v (ndarray): 1D or 2D array of y-component velocities.
     -  w (ndarray): 1D or 2D array of z-component velocities.
     -  x (ndarray): 1D or 2D array of x coordinates.
     -  y (ndarray): 1D or 2D array of y coordinates.
     -  z (ndarray): 1D or 2D array of z coordinates.
     -  method : linear, cubic, nearest
        
    Returns:
        ui, vi, wi (ndarray): velocity at interpolation points
    """
    import scipy.interpolate as scint
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()
    z = np.asarray(z).flatten()
    u = np.asarray(u).flatten()
    v = np.asarray(v).flatten()
    w = np.asarray(w).flatten()
    if algo=='griddata':
        ui = scint.griddata((x, y, z), u, (xi,yi,zi), method=method)
        vi = scint.griddata((x, y, z), v, (xi,yi,zi), method=method)
        wi = scint.griddata((x, y, z), w, (xi,yi,zi), method=method)
    elif algo=='TriInterpolator':
        tri = Triangulation(x, y, z)
        if method=='linear':
            fu = LinearTriInterpolator(tri, u)
            fv = LinearTriInterpolator(tri, v)
            fw = LinearTriInterpolator(tri, w)
        elif method=='cubic':
            fu = CubicTriInterpolator(tri, u)
            fv = CubicTriInterpolator(tri, v)
            fw = CubicTriInterpolator(tri, w)
        else:
            raise NotImplementedError()
        ui = fu(xi, yi, zi)
        vi = fv(xi, yi, zi)
        wi = fw(xi, yi, zi)
    else:
        raise NotImplementedError()
    return ui, vi, wi



if __name__ == '__main__':
    pass
