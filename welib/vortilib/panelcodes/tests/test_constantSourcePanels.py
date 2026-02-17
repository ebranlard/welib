import unittest
import numpy as np    
import os as os
import matplotlib.pyplot as plt

from welib.CFD.flows2D import flowrate2D, flowfield2D, flowfield2D_plot
from welib.vortilib.panelcodes.constantSourcePanels import CSP_solve, csp_u
from welib.vortilib.panelcodes.panel_examples import getCase

scriptDir = os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_cylinder_CpUtflowrate(self, plot=False):

        for case in ['cylinder', 'ellipse']:
            #for case in ['ellipse']:
            for method in ['CCSP']:

                # --- Case/ Geometry
                cas = getCase(case, U0=1.2)
                R = cas['R']
                
                # --- Panel method
                out = CSP_solve(cas['XP'], cas['YP'], Uxy=cas['Uxy'])

                # --- Tests
                np.testing.assert_almost_equal(out['Cp'], cas['Cp'], 3)
                np.testing.assert_almost_equal(out['Ut'], cas['Ut'], 4)

                # Flow rate should be zero
                x = R*2* np.cos(cas['theta'])
                y = R*2* np.sin(cas['theta'])
                u, v = csp_u(x, y, out['SP'], sigmas=out['sigmas'])
                Q = flowrate2D(x, y, u, v, verbose=False, ns=-1)
                np.testing.assert_almost_equal(Q, 0, 8)
                
                # --- Plots 
                if plot:
                    # --- Plot
                    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
                    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.40)
                    axes[0].plot(out['theta_CP'], out['Ut'], label='Ut')
                    axes[0].plot(out['theta_CP'], out['Un'], label='Un')
                    axes[0].plot(out['theta_CP'], cas['Ut'], 'k.', label='Theory')

                    axes[0].legend()
                    axes[1].plot(out['theta_CP'], out['Cp'], label='Cp')
                    axes[1].plot(out['theta_CP'], cas['Cp'], 'k.', label='Theory')
                    axes[1].legend()
                    
                    # --- Flow field
                    vel = lambda X, Y : csp_u(X, Y, out['SP'], out['sigmas'])
                    X, Y, U, V =  flowfield2D(vel, xmax=3.5, nx=50, U0x=cas['U0'], L=R, rel=True)
                    ax =  flowfield2D_plot(X, Y, U, V, ax=axes[2], minVal=0, maxVal=cas['maxVal'], bounded=False, rel=True)
                    ax.plot(cas['XP']/R, cas['YP']/R, 'k-',lw=3)
                    ax.set_title('Source panel, cylinder')
        if plot:
            plt.show()

if __name__ == '__main__':
    #Test().test_cylinder_CpUtflowrate(plot=True)
    unittest.main()
