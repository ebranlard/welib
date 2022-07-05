"""
References:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
#--- Legacy python 2.7
from __future__ import division
# --- General
import unittest
import numpy as np
import numpy.matlib
from scipy.special import ellipk, ellipe
# import warnings
# warnings.filterwarnings('error')


def ring_u_polar_singular(r,z,Gamma=-1,r0=1):
    """ 
    Induced velocity from a vortex ring of radius r0, located at z=0 
    Polar coordinates in and out.
    Singular formulation
    """
    # Formulation from Yoon 2004
    a = np.sqrt((r+r0)**2 + (z)**2)
    m = 4 * r * r0 / (a ** 2)
    A = (z)**2 + r**2 + r0**2
    B = - 2*r*r0
    I1 = 4.0 / a    * ellipk(m)
    I2 = 4.0 / a**3 * ellipe(m) / (1 - m)
    ur = Gamma/(4*np.pi)*r0 * z/B *(I1 - A*I2)
    uz = Gamma/(4*np.pi)*r0*( (r0 + r*A/B)*I2 - r/B*I1)
    return ur,uz

def ring_u_polar(r,z,Gamma=-1,r0=1,z0=0,epsilon=0, reg_method='Saffman'):
    """ 
    Induced velocity from a vortex ring of radius r0, located at z0, 
    Polar coordinates in, and out.
    Regularization implemented using `epsilon` and `reg_method`.  
    """
    EPSILON = 1e-07 # small value used for axis formula and singularity (when epsilon=0)
    # --- Main corpus
    r=np.asarray(r)
    z=np.asarray(z)-z0 # Note: from now on, ring is centered on 0
    ur = np.full(r.shape,np.nan)
    uz = np.full(r.shape,np.nan)

    # Enforcing  Axis formula : v_z=-Gamma/(2r0) *1 / (1+(z/r0)^2)^(3/2)  
    Iz = r < (EPSILON * r0)
    ur[Iz] = 0
    uz[Iz] = Gamma/(2*r0)*(1.0/((1 +(z[Iz]/r0)**2)**(3.0/2.0)))

    # Value on the neighborhood of the ring itself..
    if epsilon==0:
        Ir = np.logical_and(np.abs(r-r0)<(EPSILON*r0), np.abs(z)<EPSILON)
        ur[Ir]=0
        uz[Ir]=Gamma/(4*r0) # NOTE: this is arbitrary
    else:
        if reg_method=='Saffman':
            # Regularized velocity near the ring (see Saffman)
            rho = np.sqrt( (r-r0)**2 + z**2 )
            Ir = rho < epsilon/2 
            ur[Ir]=0
            uz[Ir]=Gamma/(4*np.pi*r0)*(np.log(16*r0/epsilon)-1/4) # Eq 35.36 from [1]
        elif reg_method=='linear':
            rho = np.sqrt( (r-r0)**2 + z**2 )
            Ir = rho < epsilon/2 
            # Continuity in z
            rr= r0-epsilon/2
            zz= 0
            _,uzp = ring_u_polar_singular(r0+epsilon/2, 0)
            _,uzm = ring_u_polar_singular(r0-epsilon/2, 0)
            uz_mean = uzp+uzm
            #r_ref = 
            #Ir = rho < epsilon/2 
            ur[Ir]=0
            uz[Ir]=uzp+ (uzm-uzp) * (rho[Ir]/epsilon/2)


    # --- From this point on, variables have the size of ~Iz..
    bnIz = np.logical_and(np.logical_not(Iz), np.logical_not(Ir))
    r = r[bnIz]
    z = z[bnIz]
    ur[bnIz], uz[bnIz] = ring_u_polar_singular(r,z,Gamma,r0)

    return ur,uz


def ring_u(Xcp,Ycp,Zcp,Gamma=-1,R=1,polar_out=True,epsilon=0):
    """ 
    Compute the induced velocity from a vortex ring located at z=0
    Takes cartesian coordinates as input, returns polar or cart depending on `polar_out`
    """
    EPSILON = 1e-07
    # --- Main corpus
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)
    if Xcp.shape==(0,):
        if polar_out:
            return np.array([]), np.array([])
        else:
            return np.array([]),np.array([]),np.array([]),

    r = np.sqrt(Xcp**2 + Ycp**2)
    z = Zcp
    ur = np.full(r.shape,np.nan)
    uz = np.full(r.shape,np.nan)

    # Enforcing  Axis formula : v_z=-Gamma/(2R) *1 / (1+(z/R)^2)^(3/2)  
    Iz = r < (EPSILON * R)
    ur[Iz] = 0
    uz[Iz] = Gamma/(2*R)*(1.0/((1 +(z[Iz]/R)**2)**(3.0/2.0)))

    # Value on the neighborhood of the ring itself..
    if epsilon==0:
        Ir = np.logical_and(np.abs(r-R)<(EPSILON*R), np.abs(z)<EPSILON)
        ur[Ir]=0
        uz[Ir]=Gamma/(4*R) # NOTE: this is arbitrary
    else:
        Ir = np.logical_and(np.abs(r-R)<(EPSILON*R), np.abs(z)<EPSILON)
        ur[Ir]=0
        uz[Ir]=Gamma/(4*np.pi*R)*(np.log(8*R/epsilon)-1/4) # Eq 35.36 from [1]

    # --- From this point on, variables have the size of ~Iz..
    bnIz = np.logical_and(np.logical_not(Iz), np.logical_not(Ir))
    r = r[bnIz]
    z = z[bnIz]

    # Formulation uses Formula from Yoon 2004
    a = np.sqrt((r+R)**2 + (z)**2)
    m = 4 * r * R / (a ** 2)
    A = (z)**2 + r**2 + R**2
    B = - 2*r*R
    K = ellipk(m)
    E = ellipe(m)
    I1 = np.multiply(4.0/a, K)
    I2 = 4.0 / a ** 3.0 * E / (1 - m)
    ur[bnIz] =Gamma/(4*np.pi)*R*(np.multiply(((z)/B),(I1 - np.multiply(A,I2))))
    uz[bnIz] =Gamma/(4*np.pi)*R*(np.multiply((R + np.multiply(r,A) / B),I2) - np.multiply(r/B,I1))



    # Formulation uses Formula from Yoon 2004
#     r=r[bnIz]
#     x=z[bnIz]
#     r0=R
#     x0=0
#     a  = np.sqrt((r+r0)**2 + (x-x0)**2)
#     m  = 4 * r * r0 / (a ** 2)
#     A  = (x-x0)**2 + r**2 + r0**2
#     B  = - 2*r*r0
#     I1 = 4.0/a    * ellipk(m)
#     I2 = 4.0/a**3 * ellipe(m)/(1 - m)
#     ur = Gamma/(4*np.pi)*r0*(x-x0)/B * (I1 - A*I2)
# ur[bnIz]


    if polar_out:
        return ur,uz
    else:
        psi = np.arctan2(Ycp,Xcp)     ;
        ux=ur*np.cos(psi)
        uy=ur*np.sin(psi)
        return ux,uy,uz

def rings_u(Xcp,Ycp,Zcp,Gamma_r,Rr,Xr,Yr,Zr,polar_out=True,epsilon=0):
    """ 
    Compute the induced velocity from nRings vortex rings
        nRings: number of main rings
    TODO: angles

    INPUTS: 
        Xcp,Ycp,Zcp: cartesian coordinates of control points where the velocity field is not be computed
        Gamma_t : array of size (nRings), intensity of each rings
        R       : array of size (nRings), radius of each rings
        Xr,Yr,Zr: arrays of size (nRings), center of each rings
    """
    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)
    if polar_out:
        ur = np.zeros(Xcp.shape)
        uz = np.zeros(Xcp.shape)
        for ir,(Gamma,R,xr,yr,zr) in enumerate(zip(Gamma_r,Rr,Xr,Yr,Zr)):
            if np.abs(Gamma) > 0:
                #print('Gamma',Gamma,'xr',xr,yr,zr,R)
                u1 = ring_u(Xcp-xr,Ycp-yr,Zcp-zr,Gamma,R,polar_out=polar_out,epsilon=epsilon)
                ur = ur + u1[0]
                uz = uz + u1[1]
        return ur,uz
    else:
        ux = np.zeros(Xcp.shape)
        uy = np.zeros(Xcp.shape)
        uz = np.zeros(Xcp.shape)
        for ir,(Gamma,R,xr,yr,zr) in enumerate(zip(Gamma_r,Rr,Xr,Yr,Zr)):
            if np.abs(Gamma) > 0:
                u1 = ring_u(Xcp-xr,Ycp-yr,Zcp-zr,Gamma,R,polar_out=polar_out,epsilon=epsilon)
                ux = ux + u1[0]
                uy = uy + u1[1]
                uz = uz + u1[2]
        return ux,uy,uz

# --------------------------------------------------------------------------------}
# --- TEST 
# --------------------------------------------------------------------------------{
class TestRing(unittest.TestCase):
    def test_Ring_singularities(self):
#         import warnings
#         warnings.filterwarnings('error')
        # ---- r=0, z=0, uz=Gamma/(2R)
        ur,uz=ring_u(0,0,0)
        np.testing.assert_almost_equal(uz,-1/2)
        ur,uz=ring_u_polar(0,0)
        np.testing.assert_almost_equal(uz,-1/2)
        # ---- r=1, z=0, NOTE: arbitrary
        ur,uz=ring_u(1,0,0)
        np.testing.assert_almost_equal(uz,-1/4)
        ur,uz=ring_u_polar(1,0)
        np.testing.assert_almost_equal(uz,-1/4)

    def test_Ring_axis(self):
        # Test that ring gives axis formula, eq 35.11 in reference [1]
        Gamma, R = -1, 10
        z=np.linspace(-2*R,2*R,20)
        x=z*0
        y=z*0
        ur,uz=ring_u(x,y,z,Gamma,R)
        uz_ref = Gamma/(2*R)*(1.0/(1 + (z/R)**2)**(3.0/2.0))
        np.testing.assert_almost_equal(uz,uz_ref,decimal = 7)
        np.testing.assert_almost_equal(ur,uz*0  ,decimal = 7)
        ur,uz=ring_u_polar(x,z,Gamma,R)
        np.testing.assert_almost_equal(uz,uz_ref,decimal = 7)
        np.testing.assert_almost_equal(ur,uz*0  ,decimal = 7)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(z,uz)
        #plt.plot(z,uz_ref,'k--')
        #plt.plot(z,ur)
        #plt.show()

    def test_Ring_axis_approx(self):
        # Test the approximate axis formula, close to axis, eq 35.22 in reference [1]
        Gamma, R= -1, 10
        z=np.linspace(-2*R,2*R,21)
        r=z*0+ 0.1*R
        y=z*0
        ur,uz=ring_u(r,y,z,Gamma,R)
        uz_ref =   Gamma/(2*R)*1/(1 + (z/R)**2)**(3/2)
        ur_ref = 3*Gamma/(4*R)*1/(1 + (z/R)**2)**(5/2) * r*z/R**2
        np.testing.assert_almost_equal(ur, ur_ref,decimal=3)
        np.testing.assert_almost_equal(uz, uz_ref,decimal=3)
        ur,uz=ring_u_polar(r,z,Gamma,R)
        np.testing.assert_almost_equal(ur, ur_ref,decimal=3)
        np.testing.assert_almost_equal(uz, uz_ref,decimal=3)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(z,ur_ref)
        #plt.plot(z,ur,'--')
        #plt.plot(z,uz_ref)
        #plt.plot(z,uz,'--')
        #plt.show()
        
    def test_Ring_regularized(self):
        R=1
        r=np.linspace(0,2,1000)
        z=r*0
        ur,uz = ring_u_polar(r,z, epsilon=0.1, reg_method='linear')
        zp=r*0+0.049
        zm=r*0-0.049
        urp,uzp = ring_u_polar(r,zp, epsilon=0.1, reg_method='linear')
        urm,uzm = ring_u_polar(r,zm, epsilon=0.1, reg_method='linear')
#         import matplotlib.pyplot as plt
#         plt.figure()
# #         plt.plot(r,ur,label='ur')
# #         plt.plot(r,urp,label='ur +')
# #         plt.plot(r,urm,'--',label='ur -')
# #         plt.legend
#         plt.figure()
#         plt.plot(r,uz,label='uz')
# #         plt.plot(r,uzp,label='uz +')
# #         plt.plot(r,uzm,'--',label='uz -')
#         plt.legend
#         plt.show()
        pass

    def test_Ring_doublet_approx(self):
        # Far away, the velocity field from a ring is similar to a doublet
        try:
            import welib.vortilib.elements.VortexDoublet as vd
        except:
            print('Test skipped, vortex doublet missing')
            return
        Gamma, R= -20000, 10
        mz= Gamma * np. pi * R**2 # Doublet intensity

        # --- Circular control points (centered on ring) far away from ring
        beta   = np.linspace(0,2*np.pi, 37)[1:-1]
        rProbe = 20*R
        x = rProbe*np.sin(beta)
        z = rProbe*np.cos(beta)
        y = z*0
        uxr,uyr,uzr=ring_u        (x,y,z,Gamma,R, polar_out=False)
        uxd,uyd,uzd=vd.doublet_u  (x,y,z,m=[0,0,mz])
        bNx = np.abs(uxr)>1e-3
        reldiffx=(uxr[bNx]-uxd[bNx])/(uxr[bNx])*100
        bNz = np.abs(uzr)>1e-3
        reldiffz=(uzr[bNz]-uzd[bNz])/(uzr[bNz])*100

        np.testing.assert_equal(np.max(np.abs(reldiffx))< 6, True )  
        np.testing.assert_equal(np.max(np.abs(reldiffz))< 6, True )
        np.testing.assert_almost_equal(uxr, uxd,decimal=3)
        np.testing.assert_almost_equal(uyr, uyd,decimal=8)
        np.testing.assert_almost_equal(uzr, uzd,decimal=3)

        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1,1)
        #ax.plot(beta[bNx],np.abs(reldiffx), label='reldiffx')
        #ax.plot(beta[bNz],np.abs(reldiffz), label='reldiffz')
        #ax.plot(beta,uxr    , label='x ring')
        #ax.plot(beta,uxd    , label='x doublet')
        #ax.plot(beta,uzr    , label='z ring')
        #ax.plot(beta,uzd    , label='z doublet')
        #ax.legend()
        #plt.show()
# 
    def test_Ring_rotor(self):
        pass
        # Test that induction on the rotor is constant, equal to gamma/2, see [1]
#         gamma_t, R= -1, 10
#         eps=10**-6 *R
#         x=np.linspace(0, 2*R, 1000)
#         x=np.linspace(R-eps, R+eps, 100)
#         y=x*0
#         z=x*0
#         ur,uz=ring_u(x,y,z,gamma_t,R)
        #uz_ref=[gamma_t/2]*len(x)
        #np.testing.assert_almost_equal(uz,uz_ref,decimal=7)
#         import matplotlib.pyplot as plt
#         plt.figure()
#         plt.plot(z,ur_ref)
#         plt.plot(x,ur,'--',label='ur')
#         plt.plot(z,uz_ref)
#         plt.plot(x,uz,'--',label='uz')
#         plt.legend()
#         plt.show()
# 
# 
#     def test_rings(self):
#         # Test that induction on the rotor is similat to the one from a bunch of rings
#         gamma_t, R= -1, 10
#         eps=10**-6 *R
#         x=np.linspace(-(R-eps), R-eps, 10)
#         y=x*0
#         z=x*0
#         ur,uz=vc_tang_u(x,y,z,gamma_t,R)
# #         uz_ref=[gamma_t/2]*len(x)
# #         np.testing.assert_almost_equal(uz,uz_ref,decimal=7)
#         print(ur)
#         print(uz)


if __name__ == "__main__":
#     TestRing().test_singularities()
#     TestRing().test_rotor()
    unittest.main()

