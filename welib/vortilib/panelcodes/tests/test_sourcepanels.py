import unittest
import numpy as np    
import os as os
from welib.CFD.flows2D import *
from welib.vortilib.panelcodes.sourcepanels import *

scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_flowrate(self, plot=False):
        # --- Geometry
        m  = 151
        m  = 31
        U0 = 1
        R  = 2
        theta   =-np.linspace(0,2*np.pi, m+1)
        xCP = R*np.cos(theta)
        yCP = R*np.sin(theta)
        
        # Panel method
        sigmas, out = CCSP_panel_solve(xCP, yCP, Uxy=(U0, 0))

        # --- Tests
        theta_mid = out['theta_CP']
        Ut_theory = -2.0*U0*np.sin(theta_mid)
        Cp_theory  = 1-(Ut_theory)**2/U0**2

        np.testing.assert_almost_equal(out['Cp'], Cp_theory)
        np.testing.assert_almost_equal(out['Ut'], -Ut_theory)

        # Flow rate should be zero
        x = R*2* np.cos(theta)
        y = R*2* np.sin(theta)
        u, v = ccsp_u(x, y, out['SP'], sigmas=out['sigmas'])
        Q = flowrate2D(x, y, u, v, verbose=False, ns=-1)
        np.testing.assert_almost_equal(Q, 0, 8)
        
        if plot:
            # --- Plot
            fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.40)
            axes[0].plot(out['theta_CP'], out['Ut'], label='Ut')
            axes[0].plot(out['theta_CP'], out['Un'], label='Un')
            axes[0].plot(out['theta_CP'], Ut_theory, 'k.', label='Theory')

            axes[0].legend()
            axes[1].plot(out['theta_CP'], out['Cp'], label='Cp')
            axes[1].plot(out['theta_CP'], Cp_theory, 'k.', label='Theory')
            axes[1].legend()
            
            # --- Flow field
            vel = lambda X, Y : ccsp_u(X, Y, out['SP'], out['sigmas'])
            X, Y, U, V =  flowfield2D(vel, xmax=3.5, nx=50, U0x=U0, L=R, rel=True)
            ax =  flowfield2D_plot(X, Y, U, V, ax=axes[2], minVal=0, maxVal=2, bounded=False, rel=True)
            ax.plot(xCP/R, yCP/R, 'k-',lw=3)
            ax.set_title('Source panel, cylinder')

            plt.show()

if __name__ == '__main__':
    unittest.main()
