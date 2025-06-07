import unittest
import numpy as np
import os
from welib.yams.flexibility import *
import welib.weio as weio

try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
    def test_sectionLoads_methods(self):
        # o
        #
        #
        n=20
        z = np.linspace(0,100,n)
        p = 1 + (z/100)**3
        #p=*0
        M_top=1000
        F_top=400

        F_lumped =np.zeros(p.shape)
        #F_lumped[int(n/2)]=1000
        #F_lumped[int(3*n/4)]=1000

        F_tot = trapezoid(p,z) + F_top + np.sum(F_lumped)
        M_tot = trapezoid(p*z ,z) + F_top*z[-1]

        Fs, Ms   = beamSectionLoads1D(z, p, F_top, M_top, s=1, F_lumped=F_lumped, method='plin')
        Fs2, Ms2 = beamSectionLoads1D(z, p, F_top, M_top, s=1, F_lumped=F_lumped, method='manual')
        Fs3, Ms3 = beamSectionLoads1D(z, p, F_top, M_top, s=1, F_lumped=F_lumped, method='cumtrapz')

        #import matplotlib.pyplot as plt
        #fig,axes = plt.subplots(1, 2, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        #axes[0].plot(Fs, z, label='Fs')
        #axes[0].plot(Fs2, z, '+', label='Fs (manual)')
        #axes[0].plot(Fs3, z, '.', label='Fs (cumtrapz)')
        #axes[0].plot([F_top,F_top], [z[0],z[-1]], 'k--')
        #axes[0].plot([F_tot,F_tot], [z[0],z[-1]], 'k--')

        #axes[1].plot(Ms, z, label='Ms')
        #axes[1].plot(Ms2, z, '+', label='Ms (manual)')
        #axes[1].plot(Ms3, z, '.', label='Ms (cumtrapz)')
        #axes[1].plot([M_top,M_top], [z[0],z[-1]], 'k--')
        #axes[1].plot([M_tot,M_tot], [z[0],z[-1]], 'k--')
        #
        #axes[0].set_ylabel('z [m]')
        #axes[0].set_xlabel('F [N]')
        #axes[1].set_xlabel('M [Nm]')
        #axes[0].legend()
        #axes[1].legend()
        #plt.show()

        
        dz=np.diff(z)[0]
        dMsdz = np.diff(Ms)/dz
        dFsdz = np.diff(Fs)/dz
        z2=z[0:-1]+np.diff(z)/2
        p2   = np.interp(z2, z,p)
        Fs_2 = np.interp(z2, z,Fs)

        #fig,axes = plt.subplots(1, 2, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        #axes[0].plot(p    , z ,'-', label='p')
        #axes[0].plot(p2   , z2,'-', label='p')
        #axes[0].plot(-dFsdz, z2,'--',     label='dFs/dz')

        #axes[1].plot(Fs   , z , '-', label='Fs')
        #axes[1].plot(-dMsdz, z2, '--',label='dMs/dz')
        #
        #axes[0].set_ylabel('z [m]')
        #axes[0].set_xlabel('p [N/m]')
        #axes[1].set_xlabel('F [N]')
        #plt.show()


        # Last value must equal Top
        np.testing.assert_allclose(Fs[-1] , F_top,rtol=1e-5)
        np.testing.assert_allclose(Ms[-1] , M_top,rtol=1e-5)
        # First value must equal integral
        np.testing.assert_allclose(Fs[0] , F_tot,rtol=1e-5)
        #np.testing.assert_allclose(Ms[-1] , M_tot,rtol=1e-5)
        # Methods must match
        np.testing.assert_allclose(Fs,Fs2,rtol=1e-5)
        np.testing.assert_allclose(Fs,Fs3,rtol=1e-5)
        np.testing.assert_allclose(Ms,Ms2,rtol=1e-3)
        np.testing.assert_allclose(Ms,Ms3,rtol=1e-3)

        # Section force gradient equal -p
        np.testing.assert_allclose(dFsdz, -p2  , rtol=1e-4)
        np.testing.assert_allclose(dMsdz, -Fs_2, rtol=1e-3)


if __name__=='__main__':
    #Test().test_sectionLoads_methods()
    unittest.main()
