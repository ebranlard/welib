import unittest
import numpy as np    
import os
from welib.hydro.wavekin import *

scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):
    def test_wavenumber(self):
        # Test for one period
        T = 12.  # period [s]
        g = 9.81 # gravity [m/s^2]
        h = 30.  # water depth [m]
        f = 1/T

        k= wavenumber(f, h, g)
        wavelength = 2*np.pi/k 
        np.testing.assert_almost_equal(wavelength, [177.04211]*2, 5)

        # Test for an array of periods
        T = np.array([12,12]) # period [s]
        f = 1./T

        k= wavenumber(f, h, g)
        wavelength = 2*np.pi/k 
        np.testing.assert_almost_equal(wavelength, [177.04211]*2, 5)

    def test_kinematics(self):
        g = 9.81 # gravity [m/s^2]
        h = 30.  # water depth [m]

        # --- Test for one frequency, one point, multiple time
        a    = 3    # wave peak amplitude [m]
        T    = 12.  # period [s]
        eps  = 0    # phase shift
        x, z = 0, 0 # position where kinematics are evaluated
        f    = 1./T
        k    = wavenumber(f, h, g)
        time = np.arange(0,2*T,T/101)
        vel,acc = kinematics2d(a, f, k, eps, h, time, z, x)
        eta = elevation2d(a, f, k, eps, time, x)
        np.testing.assert_almost_equal(np.max(eta),  a )
        np.testing.assert_almost_equal(np.max(vel), 2*np.pi*f*a      * np.cosh(k*(z+h)) / np.sinh(k*h), 4)
        np.testing.assert_almost_equal(np.max(acc), (2*np.pi*f)**2 *a* np.cosh(k*(z+h)) / np.sinh(k*h), 4)
        self.assertEqual(acc.shape, (len(time),) )

        # --- Test for multiple frequencies, one point, multiple time
        a    = np.array([3,  5]) # wave peak amplitude [m]
        T    = np.array([12.,9]) # period [s]
        eps  = np.array([np.pi/2, 0]) # phase shift [rad]
        x, z = 0, 0 # position where kinematics are evaluated
        f    = 1./T
        k    = wavenumber(f, h, g)
        time = np.arange(0,2*T[0],T[0]/101)
        vel,acc = kinematics2d(a, f, k, eps, h, time, z, x)
        eta = elevation2d(a, f, k, eps, time, x)
        np.testing.assert_almost_equal(np.max(eta),  7.99856, 3)
        np.testing.assert_almost_equal(np.max(vel),  5.7727, 3)
        self.assertEqual(acc.shape, (len(time),) )

        # --- Test for multiple frequencies, multiple points, one time
        a    = np.array([3,  3])  # wave peak amplitude [m]
        T    = np.array([12.,12]) # period [s]
        eps  = np.array([0 , 0])  # phase shift [rad]
        z    = np.linspace(-h, 0, 10) # position where kinematics are evaluated
        x    = z*0
        f    = 1./T
        k    = wavenumber(f, h, g)
        time = 0
        vel,acc = kinematics2d(a, f, k, eps, h, time, z, x)
        eta = elevation2d(a, f, k, eps, time, x)
        self.assertEqual(acc.shape, z.shape)

        # --- Test for multiple frequencies, multiple points, multiple times
        a    = np.array([3, ])  # wave peak amplitude [m]
        T    = np.array([12.]) # period [s]
        eps  = np.array([0  ])  # phase shift [rad]
        vz    = np.linspace(-h, 0, 10) # position where kinematics are evaluated
        vx    = np.linspace(-1,1,3)
        X,Z   = np.meshgrid(vx,vz)
        f    = 1./T
        k    = wavenumber(f, h, g)
        time = np.linspace(0,2*T[0]/2,5)
        vel,acc = kinematics2d(a, f, k, eps, h, time, Z, X)
        #eta = elevation2d(a, f, k, eps, time, x)
        np.testing.assert_array_equal(acc.shape, np.concatenate((Z.shape,time.shape)))
        np.testing.assert_almost_equal(np.max(vel[-1,1,:]), 2*np.pi*f*a      * np.cosh(k*(vz[-1]+h)) / np.sinh(k*h), 4)
        np.testing.assert_almost_equal(np.max(acc[-1,1,:]), (2*np.pi*f)**2 *a* np.cosh(k*(vz[-1]+h)) / np.sinh(k*h), 4)


    def test_wave_components(self):
        from welib.tools.stats import rsquare, mean_rel_err
        g = 9.81 # gravity [m/s^2]
        h = 30.  # water depth [m]

        # --- Test for multiple frequencies, multiple points, one time
        compFile = os.path.join(scriptDir, '../../../data/Monopile/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=50.csv')
        df = pd.read_csv(compFile)
        fp  = df['Frequency_[Hz]']
        ap  = df['Amplitude_[m]']
        epsp= df['Phase_[rad]']
        x    = 0
        kp   = wavenumber(fp, h, g)
        time = np.arange(0, 100, 0.10)
        eta_ref = elevation2d(ap, fp, kp, epsp, time, x)
        pSS = wave_components( eta_ref, time, water_depth=h, g=g, aThreshold=1e-4)
        eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], time, x)
        np.testing.assert_almost_equal(eta_sim, eta_ref, 4)


        R2,_ = rsquare(y_ref=eta_ref, y_sim=eta_sim)
        eps = mean_rel_err(y1=eta_ref, y2=eta_sim, method='meanabs')

        np.testing.assert_almost_equal(eta_sim, eta_ref, 4)

        np.testing.assert_array_less(1-R2,0.001)
        np.testing.assert_array_less(eps, 0.009)


if __name__ == '__main__':
    MyDir=os.path.dirname(__file__)
    unittest.main()
