import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


def cylinder(R1=2, R2=1, z1=0, z2=1, nz=3, nTheta=15, thetaMax=2*np.pi):
    """ 
    Generate a cylinder mesh. Optionally tapered, or incomplete (thetaMax<2*np.pi)
    With z as main axis

    Example:
        ax.plot_surface(X, Y, Z, color = 'w', rstride = 1, cstride = 1)
    """
    # Tapered Cylinder
    vtheta = np.linspace(0 ,thetaMax,nTheta)
    vz     = np.linspace(z1, z2, nz)
    vr     = np.linspace(R1, R2, nz)
    theta, z = np.meshgrid(vtheta, vz)
    _, R = np.meshgrid(vtheta, vr)
    X = R* np.cos(theta)
    Y = R* np.sin(theta)
    Z = z
    return X,Y,Z

def arbitrary_cylinder(x, y, z1=0, z2=1, nz=3, nTheta=15):
    """ 
    Generate an arbitrary (non circular) cylinder mesh
    With z as main axis
    INPUTS:
     - x,y: arrays representing the surface outline in a 2d plane

    Example:
        theta=np.linspace(0,2*np.pi,30)
        x=np.cos(theta)
        y=np.sin(2*theta)
        X,Y,Z = arbitrary_cylinder(x, y, z1=0, z2=1, nz=3, nTheta=15)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        ax.plot_surface(X, Y, Z, color = 'w', rstride = 1, cstride = 1)
    """
    assert(nz>=2)
    from welib.tools.curves import curve_interp
    # Interpolate
    x,y = curve_interp(x,y,n=nTheta)
    # Generate cylindrical surface
    vz     = np.linspace(z1, z2, nz)
    X, Z = np.meshgrid(x, vz)
    Y, _ = np.meshgrid(y, vz)
    return X,Y,Z

def torus(r=0.25, R=1, nTheta=32, nPhi=31, thetaMax=2*np.pi, phiMax=2*np.pi):
    """ 
    Generate a torus mesh. 
    For now, Torus centered about 0,0,0, with z as main axis

    Example:
        ax.plot_surface(X, Y, Z, color = 'w', rstride = 1, cstride = 1)
    """
    vtheta = np.linspace(0,thetaMax, nTheta)
    vphi   = np.linspace(0,phiMax, nPhi    )
    theta, phi = np.meshgrid(vtheta, vphi)
    X = (R + r * np.cos(phi)) * np.cos(theta)
    Y = (R + r * np.cos(phi)) * np.sin(theta)
    Z = r * np.sin(phi)
    return X,Y,Z

