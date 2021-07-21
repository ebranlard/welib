"""
Reference:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
#--- Legacy python 2.7
from __future__ import division
from __future__ import print_function
# --- General
import unittest
import numpy as np
import numpy.matlib

# --------------------------------------------------------------------------------}
# --- Semi infinite vortex line
# --------------------------------------------------------------------------------{
def vl_semiinf_u(xa,ya,za,xe,ye,ze,Gamma,visc_model,t):
    """ Induced velocity at point A, from a semi infinite vortex line starting at point 0, and directed along e
    See fUi_VortexlineSemiInfinite
    TODO vectorize 
    """
    norm_a      = np.sqrt(xa**2 + ya**2 + za**2)
    denominator = norm_a * (norm_a - (xa*xe + ya*ye + za*ze))
    crossprod   = np.array([ye*za - ze*ya, ze*xa - xe*za, xe*ya - ye*xa])
    # check for singularity */
    if (denominator < 1e-12 or norm_a < 1e-05):
        return np.array([[0],[0],[0]])
    # viscous model */
    Kv = 1.0
    if visc_model==0:
        # No vortex core model 
        Kv = 1.0
    elif visc_model==1:
        # Rankine  - t<=>rc 
        norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
        h = np.sqrt(crossprod[0]**2 + crossprod[1]**2+ crossprod[2]**2) / norm_r0
        if (h < t):
            Kv = h * h / t / t
        else:
            Kv = 1.0
    elif visc_model==2:
        # Lamb-Oseen  - t<=>rc
        norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
        h = np.sqrt(crossprod[0]**2 + crossprod[1]**2+ crossprod[2]**2) / norm_r0
        Kv = 1.0 - np.exp(- 1.25643 * h * h / t / t)
    elif visc_model==3:
        # Vatistas n=2  - t<=>rc */
        norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
        h = np.sqrt(crossprod[0]**2 + crossprod[1]**2+ crossprod[2]**2) / norm_r0
        Kv = h * h / np.sqrt(t ** 4 + h ** 4)
    elif visc_model==4:
        # Cut-off radius no dimensions - delta*norm(r0)^2 - t<=>delta^2
        norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
        denominator = denominator + t * norm_r0
        Kv = 1.0
    elif visc_model==5:
        # Cut-off radius dimension (delta l_0)^2 - t<=>(delta l)^2 */
        Kv = 1.0
        denominator = denominator + t
    Kv = Gamma*Kv /4.0/ np.pi / denominator
    Uout = Kv * crossprod
    return Uout

def vl_semiinf_straight_u_polar(vr,vpsi,vz,Gamma_r=-1,polar_out=False):
    """
    Induced velocity from a semi-infinite line along the z axis (e.g. root vortex)
    Takes polar coordinates as inputs, returns velocity either in Cartesian (default) or polar.
    INPUTS:
       vr,vpsi,vz : control points in polar coordinates, may be of any shape
       Gamma_r    : Root vortex circulation, negative for a wind turbine
       m =tan(chi): tangent of wake skew angle
    Reference: [1,2]"""
    EPSILON_AXIS=1e-7; # relative threshold for using axis formula
    if Gamma_r==0:
        return vr*0,vr*0,vr*0
    # Flattening
    vr       = np.asarray(vr)
    shape_in = vr.shape
    vpsi     = np.asarray(vpsi)
    vz       = np.asarray(vz)

    u_psi = Gamma_r/(4*np.pi*vr) * ( 1 + vz/np.sqrt(vr**2 + vz**2) ) 

    u_z = np.zeros(shape_in)
    if polar_out:
        u_r = np.zeros(shape_in)
        u_psi = u_psi.reshape(shape_in) 
        return (u_r,u_psi,u_z)
    else:
        u_x   = -np.sin(vpsi)*u_psi
        u_y   =  np.cos(vpsi)*u_psi
        return (u_x,u_y,u_z)

def vl_semiinf_straight_u(Xcp,Ycp,Zcp,Gamma_r=-1,polar_out=False):
    """ Like vl_semiinf_straight_u_polar, with Cartesian inputs"""
    vr, vpsi = np.sqrt(Xcp**2+Ycp**2), np.arctan2(Ycp,Xcp) # polar coords
    return  vl_semiinf_straight_u_polar(vr,vpsi,Zcp,Gamma_r,polar_out) # ux,uy,uz OR ur,upsi,uz

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestVortexLine(unittest.TestCase):
    def test_VL_semi_inf(self):
        import warnings
        warnings.filterwarnings('error')

        # --- Semi-infinite line (along a vector e) and straight analytical one should return the same
        e = np.array([0,0,1])
        r=1
        psi=np.pi/6
        z=-10
        Gamma=1
        u = vl_semiinf_u(r*np.cos(psi),r*np.sin(psi),z,e[0],e[1],e[2],Gamma,visc_model=0,t=0)
        us = vl_semiinf_straight_u_polar(r,psi,z,Gamma,polar_out=False)
        np.testing.assert_almost_equal(u,us,decimal=16)


        # Typical case for a root vortex
        r=1
        psi=0
        z=10
        Gamma=-1
        u = vl_semiinf_u(r*np.cos(psi),r*np.sin(psi),z,e[0],e[1],e[2],Gamma,visc_model=0,t=0)
        us = vl_semiinf_straight_u_polar(r,psi,z,Gamma,polar_out=False)
        np.testing.assert_equal(u[1]<0,True)
        np.testing.assert_equal(us[1]<0,True)
        #print(u)
        #print(us)


if __name__ == "__main__":
    unittest.main()


