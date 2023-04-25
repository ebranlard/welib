import unittest
import numpy as np    
import os as os
from welib.airfoils.section import Section

scriptDir = os.path.dirname(__file__)

class TestSection(unittest.TestCase):
    def test_section_x(self):
        # Structure
        m = 1
        c = 0
        k = 61.7
        x0   = -0.02 # initial position
        xd0  = 0     # initial velocity

        # Inputs
        Uy=2    # Incoming wind speed "Omega r" [m/s]
        Ux=0    # Incoming wind speed "U0" [m/s]
        theta_p = -15


        # TODO Interface likely to change
        sec = Section()
        sec._chord = 0.2
        sec._rho = 1.225
        sec._tau = 0.08
        sec._y_AQ = 0
        sec._y_AT = 0
        M = sec.setMassMatrix(m=m)
        C = sec.setDampMatrix(cxx=c, cyy=c)
        K = sec.setStifMatrix(kxx=k, kyy=k)
        sec.polarFromCSV(polarFilename=os.path.join(scriptDir,'../data/tjaere11_ds.csv'))
        sec.setParameters(sx_sim='x', ds_model='oye', di_model=None)
        sec.setConstantInputs(Ux, Uy, theta_p=theta_p*np.pi/180)

        # --- Time integration
        fs_i= 0.88792530768895 # TODO
        t = np.linspace(0, 2, 7*60)
        q0 = sec.setInitialConditions(x=x0, fs=fs_i)
        res, df = sec.integrate(t, y0=q0)
        np.testing.assert_almost_equal(df['x1'].values[-1], 0.035058, 5)
        np.testing.assert_almost_equal(df['x2'].values[-1], 0.069730, 5)
        np.testing.assert_almost_equal(df['x3'].values[-1], 0.779248, 5)
        #sec.plot_states()
        #import matplotlib.pyplot as plt
        #plt.show()

        # --- Frequencies
        u0 = [Ux, Uy, theta_p*np.pi/180]
        xop = sec.equilibrium(x0=q0, u0=u0)
        freq_d, zeta, Q, freq_0 = sec.eigA(x0=xop, u0=u0)
        #print('Frequencies',freq_0)
        #print('Damping    ',zeta*100)
        # print('Q\n',np.abs(Q))
        np.testing.assert_almost_equal(freq_0[0]  , 1.207931, 5)
        np.testing.assert_almost_equal(zeta[0]*100, 1.080974, 5)


if __name__ == '__main__':
    unittest.main()
