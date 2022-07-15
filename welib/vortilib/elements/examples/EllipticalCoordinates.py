from welib.vortilib.elements.SourceEllipsoid import *
import numpy as np
import matplotlib.pyplot as plt

def main():
    # --- Ellipse parameters
    a = 3  # semimajor axis
    b = 2  # semiminor axis 
    e = np.sqrt(a**2 - b**2)/a # eccentricity
    # e = 1/a # Values such that a*e=1
    # b= a * np.sqrt(1-e**2)

    k=a*e
    print('a    ',a)
    print('b    ',b)
    print('e    ',e)
    print('k    ',k)
    print('zeta0',1/e)

    # --- Coordinates of the ellipse
    # In elliptical coordinates: zeta=zeta0, mu in [-1,1]
    ne = 35
    zeta0 = np.ones(ne)*1/e      
    mu0   = np.linspace(-1,1,ne)
    # In Cartesian Coordinates:
    x0,y0,z0 = T_semielliptic_cart(mu0,zeta0,mu0*0,k) 
    # Cartesian coordinates, directly
    x0e = np.linspace(-a,a,ne)
    y0e = b*np.sqrt(1-(x0e/a)**2)
    # Transform back to elliptical coordinates
    #mu02,zeta02,theta02 = T_cart_semielliptic(x0,y0,z0,k)

    # --- Grid in elliptical coordinates
    nmu = 21
    nzeta=11
    vMu   = np.cos(np.linspace(-np.pi,np.pi,nmu))   # full range [-1,1]
    vZeta = np.cosh(np.linspace(0,2,nzeta))
    # vMu   = np.linspace(-1,1,nmu)   # full range [-1,1]
    # vZeta = np.unique(np.concatenate((np.linspace(1,zeta0,nzeta),np.linspace(zeta0,2*zeta0,nzeta))))
    # vZeta = np.array([zeta0*1.2])
    # Mesh grid
    MU,ZETA = np.meshgrid(vMu, vZeta)
    THETA = MU*0
    # Transform the grid in Cartesian coordinates for plotting
    X,Y,Z    = T_semielliptic_cart(MU,ZETA,THETA,k)
    X1,Y1,Z1 = T_semielliptic_cart(MU,ZETA,THETA+np.pi,k) # Add pi for lower domain (y<0)
    # Transform back to elliptcial coordinates
    MU2,ZETA2,THETA2 = T_cart_semielliptic(X,Y,Z,k)
    np.testing.assert_almost_equal(np.min(np.abs(MU2-MU)),0)
    np.testing.assert_almost_equal(np.max(np.abs(MU2-MU)),0)
    np.testing.assert_almost_equal(np.min(np.abs(ZETA2-ZETA)),0)
    np.testing.assert_almost_equal(np.max(np.abs(ZETA2-ZETA)),0)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, hspace=0.20, wspace=0.20)
    ax.plot(X    , Y   , '-', c=(0.5,0.5,0.5), lw=0.5) # "radial" grid lines    (mu=cst)
    ax.plot(X.T  , Y.T , '-', c=(0.5,0.5,0.5), lw=0.5) # "azimuthal" grid lines (zeta=cst)
    ax.plot(X1   , Y1  , '-', c=(0.5,0.5,0.5), lw=0.5) # radial, lower domain
    ax.plot(X1.T , Y1.T, '-', c=(0.5,0.5,0.5), lw=0.5) # azimuthal, lower domain
    ax.plot(x0e,y0e ,'k-',lw=2)
    ax.plot(x0e,-y0e,'k-',lw=2)
    #plt.plot(x0,y0   ,'b:',lw=2)
    ax.plot(a,0,'kd')
    ax.plot(0,b,'kd')
    ax.plot(-e*a,0,'ko')
    ax.plot( e*a,0,'ko')
    ax.axis('equal')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Elliptical Coordinates')


if __name__ == '__main__':
    main()
    plt.show()
if __name__=="__test__":
    pass
if __name__=="__export__":
    main()
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)


