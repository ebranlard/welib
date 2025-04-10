"""
Trailed (longitudinal) Vorticity from a skewed vortex cylinder.
This test reproduces the plots of the section 3.3 of the following reference:
    [1] E. Branlard, M. Gaunaa - Cylindrical vortex wake model: skewed cylinder, application to yawed or tilted rotors - Wind Energy, 2015


Coordinate system:
 rotor plane x-y  (y is vertical, x longi)
 wake is z (if not yawed)
""" 

## Parameters
# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
# --- Local
from welib.vortilib.elements.VortexCylinderSkewed import fV_Root
from welib.tools.colors import *

def pol2cart(rho,phi):
    x= rho * np.cos(phi)
    y= rho * np.sin(phi)
    return (x,y)

def main(test=False):
    #np.warnings.filterwarnings('error')
    ## Params
    Gamma_tot = 1
    Gamma_r = - 1
    visc_model = 0
    t = 0
    R = 1
    bComputeGrad = 0

    ## --- n lines along a circle - AZIMUTHAL SURVEY - WITH ROOT VORTEX
    chi    = 30 * pi / 180
    e      = np.array([np.sin(chi),0,np.cos(chi)])
    n      = 100
    vtheta = np.linspace(0,2 * pi,n)
    npsi   = 100
    vpsi   = np.linspace(0,2 * pi,npsi)
    vr     = np.array([0.2,0.5,1])
    nr     = len(vr)
    vz0    = 0
    m      = np.tan(chi)
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    v_z_root ,v_psi_root  = fV_Root(VR,VPSI,VZ0,m,Gamma_r,nout = 2)
    v_z_root0,v_psi_root0 = fV_Root(VR,VPSI,VZ0,0,Gamma_r,nout = 2)

    plt.figure()
    for ir,r in enumerate(vr):
        plt.plot(vpsi*180/pi,np.squeeze(v_z_root[:,ir,:])*4*pi*R/Gamma_tot, label='$r/R=%2.1f$'%r)
    plt.plot(vpsi*180/pi, np.squeeze(v_z_root0[:,0,:])*4*pi*R/Gamma_tot,'k--',linewidth=2,label=r'Value for $\chi=0$')
    plt.legend()
    plt.ylabel(r'$u_{z,r} 4 \pi R / \Gamma_{tot}$')
    plt.xlabel(r'Azimuthal position $\psi$ [deg]')
    plt.title(r'RootVortexOnly Axial velocity for various radii Chi%2d'%np.round(chi * 180 / pi))
    plt.xlim(np.array([0,360]))

    plt.figure()
    for ir,r in enumerate(vr):
        plt.plot(vpsi*180/pi,np.squeeze(v_psi_root[:,ir,:])*4*pi*R/Gamma_tot, label='$r/R=%2.1f$'%r)
    plt.plot(vpsi*180/pi, np.squeeze(v_psi_root0[:,0,:])*4*pi*R/Gamma_tot,'k--',linewidth=2,label=r'Value for $\chi=0$')
    plt.legend()
    plt.ylabel(r'$u_{\psi,r} 4 \pi R / \Gamma_{tot}$')
    plt.xlabel(r'Azimuthal position \psi [deg]')
    plt.title(r'RootVortexOnlyTangential velocity for various radii Chi%2d'%np.round(chi*180/pi))
    plt.ylim(np.array([-8.8,0]))
    plt.xlim(np.array([0,360]))

    ## --- Contour plot of axial velocity
    n_azimuth = 180
    n_radial  = 25
    vpsi      = np.linspace(0,360,n_azimuth) * pi / 180
    vr        = np.linspace(0.0,1.5,n_radial) * R
    vz0       = 0
    chi       = 30 * pi / 180
    m         = np.tan(chi)
    VR,VPSI,VZ0 = np.meshgrid(vr,vpsi,vz0)
    VR   = VR.squeeze()
    VPSI = VPSI.squeeze()
    VZ0  = VZ0.squeeze()
    Xplane,Yplane = pol2cart(VR,VPSI)
    # Ui
    uz = fV_Root(VR,VPSI,VZ0,m,Gamma_r,nout=1)
    # Plot
    plt.figure()
    ContourVal = np.array([0,0.01,0.05,0.1,0.2])
    ContourVal = np.unique(np.sort(np.concatenate((-ContourVal,ContourVal))/ Gamma_tot))
    CS = plt.contour(Xplane/R, Yplane/R, uz/Gamma_tot,levels=ContourVal)
    manual_locations=[ (0.1667575,0.418641), (0.1492567,- 0.4228), (0.0758036,- 0.2096), (0.0213203,0.203710), (- 0.12992,0.702902), (- 0.13028,- 0.7035), (- 0.66038,0.175041), (- 0.66048,- 0.1747)]
    lbs=plt.clabel(CS,inline=1,fontsize=10,manual=manual_locations)
    rots=[ - 350.3, - 6.172, - 6.064, - 12.97, - 36.47, - 323.4, - 336.0, - 23.83]
    for l,r in zip(lbs,rots):
        l.set_rotation(r)
    plt.axis('equal')
    vtheta = np.linspace(0,2 * pi,100)
    plt.plot(R * np.cos(vtheta),R * np.sin(vtheta),'k')
    plt.axis('square')
    plt.xlim(np.array([  1.3,-1.3])) # NOTE: flipped
    plt.ylim(np.array([- 1.3,1.3]))
    plt.title('RootVortexOnlyContourAxial_chi%03d'%np.round(np.arctan(m) * 180 / pi))
    plt.xlabel('$x/R$')
    plt.ylabel('$y/R$')

    if not test:
        plt.show()
    else:
        plt.close('all')

class Test(unittest.TestCase):
    def test_Article_Skew_RootVortex(self):
        import sys
        if sys.version_info >= (3, 0):
            main(test=True)
        else:
            print('Test skipped due to travis display error')

if __name__ == "__main__":
    main()
