"""
References:
    [1] Lamb - Hydrodynamics, p.141
    [2] Durand - Aerodynamic theory, Vol.1,  p.278
    [3] Anderson, Branlard, Vijayakumar, Johnson - Investigation of the nacelle blockage effect for a downwind turbine, Journal of Physics


NOTES on semi-elliptic coordinates (mu, zeta, omega) :
    mu, zeta: hyperboid and ellipdoid of revolution coordinates (confocal with current ellipse)
    a, e : large axis and eccentricity of the ellipse
    x     = k * mu zeta
    r     = k * sqrt(1-mu^2) sqrt(zeta^2-1)
    omega 
    e     = sqrt(a^2 + b^2)/a
    zeta0 = 1/e

    x= k \cos\theta \cosh eta = k * mu * zeta
    r= k \sin\theta \sinh eta = k * sqrt(1-mu**2) * sqrt(zeta**2 -1)

"""

# --- General
import unittest
import numpy as np
# --- Local


def T_semielliptic_cart(mu,zeta,theta,k=1):
    """ 
    Tranforms semi-elliptic coordinates to cartesian coordinates (see [3])
        k=a*e , grand axe times eccentricity
    """
    x = k*mu*zeta
    r = k*np.sqrt(1-mu**2)*np.sqrt(zeta**2-1)
    y = r*np.cos(theta)
    z = r*np.sin(theta)
    return x,y,z

def T_cart_semielliptic(x,y,z,k=1):
    """ 
    Tranforms cartesian coordinates to semi-elliptic coordinates (see [3])
        k=a*e , grand axe times eccentricity
    """
    r               = np.sqrt(y**2+z**2)
    theta           = np.arctan2(z,y)
    a               = k**2
    b               = -(k**2+x**2+r**2)
    c               = x**2
    delta           = b**2 - 4 * a*c
    zeta2           = (- b + np.sqrt(delta)) / (2*a)
    zeta            = np.sqrt(zeta2)
    mu = x/(k*zeta)
    zeta[zeta<= 1] = 1
    mu  [mu>= 1 ]  = 1
    mu  [mu<= -1]  = -1
    return mu,zeta,theta

def ellipse_coord(a,b,ne=100):
    theta = np.linspace(-np.pi,np.pi,ne)
    xe    = a * np.cos(theta)
    ye    = b * np.sin(theta)
    return xe, ye


def ser_phi_elliptic(mu,zeta,U0,a,b):
    """ 
    Potential (phi) from an Source Ellipsoid of Revolution (SER) with inputs expressed using
    ovary semi elliptic coordinates.  (see documentation at top of module)
    INPUTS: 
      - mu, zeta: hyperboid and ellipdoid of revolution coordinates (confocal with current ellipse)
      - a, e : large axis and eccentricity of the ellipse
      - U : free stream / velocity of the ellipsoid

    OUTPUTS:
      - phi: velocity potential in this semi-elliptci coordinates
    NOTE: no freestream included! 
    """
    e   = np.sqrt(a**2 - b**2)/a
    A   = U0 * a /( 1/(1-e**2) - 1/(2*e) * np.log((1+e)/(1-e)) )
    phi =np.zeros(mu.shape)
    bOK  = zeta>1
    mu   = mu[bOK]
    zeta = zeta[bOK]
    phi[bOK] = A * mu * (1/2 * zeta * np.log( (zeta+1)/(zeta-1) )  -1 )
    # phi = A * mu * (1/2 * zeta * np.log( (zeta+1)/(zeta-1) )  -1 )
    return phi

def ser_dphi_elliptic(mu,zeta,U0,a,b):
    """ Returns gradient of potential of SER in semi-elliptic coordinates
    """ 
    e   = np.sqrt(a**2 - b**2)/a
    A   = U0 * a /( 1/(1-e**2) - 1/(2*e) * np.log((1+e)/(1-e)) )
    dphi_dmu   = np.zeros(mu.shape)
    dphi_dzeta = np.zeros(mu.shape)
    bOK  = zeta>1
    mu   = mu[bOK]
    zeta = zeta[bOK]
    LOGZ = np.log((zeta+1)/(zeta-1))

    dphi_dmu[bOK]   = A *(1/2 * zeta * LOGZ  -1 )
    dphi_dzeta[bOK] = A * mu * 1/2 * ( LOGZ  - 2*zeta/(zeta**2-1) ) 
    return dphi_dmu,dphi_dzeta


def ser_psi_elliptic(mu,zeta,U0,a,b):
    """ 
    Stream function for a Source Ellipoid of Revolution (SER), see ser_phi_elliptic for documentation.
    NOTE: no freestream included! 
    """
    e   = np.sqrt(a**2 - b**2)/a
    A   = U0 * a /( 1/(1-e**2) - 1/(2*e) * np.log((1+e)/(1-e)) )
    k   = a*e
    psi =np.zeros(mu.shape)
    bOK  = zeta>1
    mu   = mu[bOK]
    zeta = zeta[bOK]
    psi[bOK] = 1/2 * A * k * (1-mu**2) * (zeta**2-1) * (1/2 * np.log( (zeta+1)/(zeta-1)) - zeta/(zeta**2 -1) ) 
    return psi

def ser_phi(X,Y,U0,a,b):
    """ Velocity potential for Source Ellipsoid of Revolution (see ser_phi_elliptic) 
    NOTE: no freestream included! 
    """
    MU,ZETA,THETA = T_cart_semielliptic(X,Y,X*0, a*np.sqrt(a**2-b**2)/a)
    phi = ser_phi_elliptic(MU,ZETA,U0,a,b)
    return phi

def ser_psi(X,Y,U0,a,b):
    """ Streamfunction for Source Ellipsoid of Revolution (see ser_psi_elliptic) 
    NOTE: no freestream included!
    """
    MU,ZETA,THETA = T_cart_semielliptic(X,Y,X*0,a*np.sqrt(a**2-b**2)/a)
    psi = ser_psi_elliptic(MU,ZETA,U0,a,b)
    return psi

def ser_u(x,y,U0,a,b):
    """ Velocity induced by Source Ellipsoid of Revolution (see ser_psi_elliptic) 
    NOTE: no freestream included! Add U0 to U at the end.
    TODO: true cartesian coordinates 
    """
    k = a*np.sqrt(a**2-b**2)/a
    mu,zeta,theta = T_cart_semielliptic(x,y,x*0,a*np.sqrt(a**2-b**2)/a)
    dphi_dmu,dphi_dzeta =  ser_dphi_elliptic(mu,zeta,U0,a,b)
    r = np.abs(y)

    k2     = k**2
    x2     = x**2
    B      = (-k2 - r**2 - x2)
    SQ_DELm = (B**2 -4*k2*x2 )**(-0.5) # inverse of square root of determinant
    XKZ2   = -x/(k*zeta**2)

    dzeta_dx =      x*SQ_DELm * (zeta-1/zeta)
    dzeta_dr = r*zeta*SQ_DELm
    dmu_dx   = XKZ2 * dzeta_dx + 1/(k*zeta)
    dmu_dr   = XKZ2 * dzeta_dr
#     AA     = (-B + SQ_DEL)
#     AAm12   = AA**(-0.5)
#     T1     = x* np.sqrt(2)*(AA/SQ_DEL -2*k2/SQ_DEL) * AAm12
#     T2     = r* np.sqrt(2)*(AA/SQ_DEL             ) * AAm12
#     dzeta_dx =     T1/(2*k)
#     dzeta_dr =     T2/(2*k)
#     dmu_dx   = -x *T1/AA + np.sqrt(2)*AAm12
#     dmu_dr   = -x *T2/AA

    u = dphi_dmu * dmu_dx + dphi_dzeta * dzeta_dx  # dphi_dx
    v = dphi_dmu * dmu_dr + dphi_dzeta * dzeta_dr  # dphi_dr

    v[y<0]*=-1

    return u,v

def ser_u_numerical(X,Y,vx,vy,U0,a,b):
    """ Velocity induced by Source Ellipsoid of Revolution (see ser_psi_elliptic) 
    NOTE: no freestream included! Add U0 to U at the end.
    TODO: true cartesian coordinates 
    """
    MU,ZETA,THETA = T_cart_semielliptic(X,Y,X*0,a*np.sqrt(a**2-b**2)/a)
    PHI = ser_phi_elliptic(MU,ZETA,U0,a,b)
    Ux=np.gradient(PHI,vx,axis=1)
    Ur=np.gradient(PHI,vy,axis=0) # NOTE: not true cartesian coordinates
    return Ux,Ur



# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestSourceEllipsoid(unittest.TestCase):
    def test_SER_Transformations(self):
        import warnings
        warnings.filterwarnings('error')
        # --- Ellipse parameters
        a = 1
        b = 0.5*a
        e = np.sqrt(a**2 - b**2)/a
        k = a*e # a * e
        zeta0 =1/e
        # --- From elliptic to cartesian and back
        nmu   = 7
        nzeta = 5
        vMu   = np.cos(np.linspace(-np.pi,np.pi,nmu))   # full range [-1,1]
        vZeta = np.cosh(np.linspace(0,2,nzeta))
        MU,ZETA = np.meshgrid(vMu, vZeta)
        THETA = MU*0 # TODO theta none 0
        X,Y,Z    = T_semielliptic_cart(MU,ZETA,THETA,k)
        MU2,ZETA2,THETA2 = T_cart_semielliptic(X,Y,Z,k)
        np.testing.assert_almost_equal(MU  ,MU2  ,10)
        np.testing.assert_almost_equal(ZETA,ZETA2,10)
        np.testing.assert_almost_equal(THETA,THETA2,10)

        # --- From cartesian to elliptic and back
        nx  = 7
        ny  = 5
        vx  = np.linspace(-2*a,2*a,nx)
        vy  = np.linspace(-2*b   ,2*b,ny)
        X,Y = np.meshgrid(vx, vy)
        MU,ZETA,THETA = T_cart_semielliptic(X,Y,Z,k)
        X2,Y2,Z2    = T_semielliptic_cart(MU,ZETA,THETA,k)
        np.testing.assert_almost_equal(X,X2,10)
        np.testing.assert_almost_equal(Y,Y2,10)
        np.testing.assert_almost_equal(Z,Z2,10)

    def test_SER_NumericalDiff(self):
        # --- 
        nx=100
        ny=nx+1
        # Ellipse parameters
        a  = 1
        b  = 0.5*a
        U0 = 10
        # --- Velocity field on grid
        vx = np.linspace(-2*a,2*a,nx)
        vy = np.linspace(-2*b   ,2*b,ny)
        X,Y = np.meshgrid(vx, vy)
        U1,V1 = ser_u_numerical(X,Y,vx,vy,U0, a, b)
        U ,V   = ser_u         (X,Y,      U0, a, b )
        bInEllipse=(X**2/a**2+Y**2/b**2)<1.01
        U[bInEllipse]  = 0
        U1[bInEllipse] = 0
        V[bInEllipse]  = 0
        V1[bInEllipse] = 0
        np.testing.assert_almost_equal(U/U0,U1/U0,1)
        np.testing.assert_almost_equal(V/U0,V1/U0,1)
        # errmax=vx[1]-vx[0]
        # fig,ax = plt.subplots(1,2)
        # im0=ax[0].pcolormesh(X, Y, np.abs(U-U1)/U0, vmin=0,vmax=errmax)
        # ax[0].plot(xe,ye,'k',lw=3)
        # fig.colorbar(im0,ax=ax[0])
        # im1=ax[1].pcolormesh(X, Y, np.abs(V-V1)/U0, vmin=0,vmax=errmax)
        # ax[1].plot(xe,ye,'k',lw=3)
        # fig.colorbar(im1,ax=ax[1])
        # plt.show()



if __name__ == "__main__":
    #     TestSourceEllipsoid().test_SER_Transformations()
    unittest.main()
